%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run LittleDog across a tilted terrain using inverse kinematics planner;

display('Clearing workspace')
clc; clear all; close all;

display('checking for dependencies')
checkDependency('iris');
checkDependency('lcmgl');
path_handle = addpathTemporary(fileparts(mfilename('fullpath')));

display('Loading robot model')
robot_options = struct();
robot_options = applyDefaults(robot_options, struct('use_bullet', true, 'terrain', RigidBodyFlatTerrain,...
                                                    'floating', true, 'ignore_self_collisions', true,...
                                                    'enable_fastqp', false, 'ignore_friction', true, 'dt', 0.001));
% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

% construct robot model
r = LittleDog('../../../LittleDog/src/LittleDog/LittleDog.urdf',robot_options);
r = r.removeCollisionGroupsExcept({'feet'});
r = compile(r);

% set initial
% state to fixed point
xstar = r.getNomState();
xstar = r.resolveConstraints(xstar);

r = r.setInitialState(xstar);
x0 = xstar;

x0.base_x = xstar.base_x + 0.85;
x0.base_z = xstar.base_z + 0.5;

xf = x0;

xf.base_x = xstar.base_x + 2;
xf.base_z = xstar.base_z + 0.5;
xf.base_yaw = x0.base_yaw + pi/8;
nq = r.getNumPositions();

display('loading obstacles');

box_size = [0.6, 2, .5];

box_tops = [0.85, 0, .5, 0,0;
            2.15, 0, .5, 0,0]';

safe_regions = iris.TerrainRegion.empty();

for j = 1:size(box_tops, 2)
  rpy = [box_tops(5,j);0;box_tops(4,j)];
  b = RigidBodyBox(box_size, box_tops(1:3,j) + [0;0;-box_size(3)/2], rpy);
  r = r.addGeometryToBody('world', b);
  offset = rpy2rotmat(rpy) * [0;0;box_size(3)/2];
  width = 0.9*box_size(1)/2;
  len = 0.9*box_size(2)/2;
  [A, b] = poly2lincon(box_tops(1,j) + offset(1) + [-width, width, width, -width],...
                       box_tops(2,j) + offset(2) + [-len, -len, len, len]);
  [A, b] = convert_to_cspace(A, b);
  normal = rpy2rotmat(rpy) * [0;0;1];
  safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], box_tops(1:3,j), normal);
end

box_size = [0.25, 0.15, .5];

box_tops = [1.3, -0.3, .5, 0,0;
            1.3, 0.3, .5, 0,0;
            1.3, -0.1, .5, 0,0;
            1.3, 0.1, .5, 0,0;
            1.5, -0.2, .5, 0,0;
            1.5, 0.2, .5, 0,0;
            1.5, 0, .5, 0,0;
            1.7, -0.1, .5, 0,0;
            1.7, 0.1, .5, 0,0;
            1.7, 0.3, .5, 0,0;
            1.7, -0.3, .5, 0,0]';

for j = 1:size(box_tops, 2)
  rpy = [box_tops(5,j);0;box_tops(4,j)];
  b = RigidBodyBox(box_size, box_tops(1:3,j) + [0;0;-box_size(3)/2], rpy);
  r = r.addGeometryToBody('world', b);
  offset = rpy2rotmat(rpy) * [0;0;box_size(3)/2];
  width = 0.85*box_size(1)/2;
  len = 0.85*box_size(2)/2;
  [A, b] = poly2lincon(box_tops(1,j) + offset(1) + [-width, width, width, -width],...
                       box_tops(2,j) + offset(2) + [-len, -len, len, len]);
  [A, b] = convert_to_cspace(A, b);
  normal = rpy2rotmat(rpy) * [0;0;1];
  safe_regions(end+1) = iris.TerrainRegion(A, b, [], [], box_tops(1:3,j), normal);
end

r = compile(r);

display('setting terrain');
height_map = RigidBodyHeightMapTerrain.constructHeightMapFromRaycast(r,x0(1:nq),-1:.015:4,-1:.015:1,10);
r = r.setTerrain(height_map).compile();
r.default_walking_params.drake_min_hold_time = 1.0;

display('Building visualizer');
v = r.constructVisualizer();
v.display_dt = 0.01;
v.draw(0,x0(1:nq))

display('Drawing goal');
lcmgl1 = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_pose');
lcmgl1.glColor3f(0,0,0);
f_height = xf.base_z - xstar.base_z;
lcmgl1.glPushMatrix();
lcmgl1.circle(xf.base_x, xf.base_y, f_height, 0.15);
lcmgl1.glPopMatrix();
lcmgl1.switchBuffers();

lcmgl3 = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_rot');
lcmgl3.glColor3f(0,0,0);
lcmgl3.glTranslated(xf.base_x, xf.base_y, f_height);
axis = rpy2axis([0,0,xf.base_yaw]);
axis = axis([4,1:3]); 
axis(1) = 180 / pi * axis(1);
lcmgl3.glRotated(axis(1), axis(2), axis(3), axis(4));
lcmgl3.glPushMatrix();
len = 0.15;
lcmgl3.glTranslated(len / 2, 0, 0);
lcmgl3.drawArrow3d(len, 0.02, 0.02, 0.005);
lcmgl3.glPopMatrix();
lcmgl3.switchBuffers();

display('planning footsteps')
footstep_plan = r.planFootstepsAdditive(x0(1:nq), xf(1:nq), safe_regions, 4, 12);
lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'footstep_plan');
  footstep_plan.draw_lcmgl(lcmgl);
  lcmgl.switchBuffers()
display('fitting footsteps to terrain')
nsteps = length(footstep_plan.footsteps);
for j = 7:nsteps
  if ~footstep_plan.footsteps(j).pos_fixed(3)
    footstep_plan.footsteps(j) = Quad_fitStepToTerrain(r, footstep_plan.footsteps(j));
  end
end

display('Adding terrain profiles')
for j = 7:length(footstep_plan.footsteps)
  [~, contact_width] = Quad_contactVolume(r, footstep_plan.footsteps(j-4), footstep_plan.footsteps(j));
  footstep_plan.footsteps(j).terrain_pts = Quad_sampleSwingTerrain(r, footstep_plan.footsteps(j-4), footstep_plan.footsteps(j), contact_width/2, struct());
end

display('Planning ZMP trajectory')
settings = Quad_QPLocomotionPlanSettings.fromQuadrupedFootstepPlan(footstep_plan, r, x0(1:nq));

%display('generating walking plan')

display('Drawing plan')
if isa(v, 'BotVisualizer')
  lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'footstep_plan');
  footstep_plan.draw_lcmgl(lcmgl);
  lcmgl.switchBuffers();
  lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'walking_plan');
  settings.draw_lcmgl(lcmgl);
  lcmgl.switchBuffers();
else
  figure(25)
  footstep_plan.draw_2d();
end

display('Running Inverse Kinematics')
[xtraj, htraj, ts] = Quad_planWalkingStateTraj(r, settings);


display('Running plan')
v.playback(xtraj, struct('slider', true));

display('Drawing goal');
lcmgl1 = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_pose');
lcmgl1.glColor3f(0,0,0);
f_height = xf.base_z - xstar.base_z;
lcmgl1.glPushMatrix();
lcmgl1.circle(xf.base_x, xf.base_y, f_height, 0.15);
lcmgl1.glPopMatrix();
lcmgl1.switchBuffers();

lcmgl3 = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_rot');
lcmgl3.glColor3f(0,0,0);
lcmgl3.glTranslated(xf.base_x, xf.base_y, f_height);
axis = rpy2axis([0,0,xf.base_yaw]);
axis = axis([4,1:3]); 
axis(1) = 180 / pi * axis(1);
lcmgl3.glRotated(axis(1), axis(2), axis(3), axis(4));
lcmgl3.glPushMatrix();
len = 0.15;
lcmgl3.glTranslated(len / 2, 0, 0);
lcmgl3.drawArrow3d(len, 0.02, 0.02, 0.005);
lcmgl3.glPopMatrix();
lcmgl3.switchBuffers();