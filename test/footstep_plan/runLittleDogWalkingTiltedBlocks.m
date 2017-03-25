%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run LittleDog across a tilted terrain using inverse kinematics planner

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
r = LittleDog('../../src/LittleDog/LittleDog.urdf',robot_options);
r = r.removeCollisionGroupsExcept({'feet'});
r = compile(r);

% set initial
% state to fixed point
xstar = r.getNomState();
xstar = r.resolveConstraints(xstar);

r = r.setInitialState(xstar);
x0 = xstar;
x0.base_z = xstar.base_z + 0.03;
x0.base_pitch = x0.base_pitch + 3*pi/180;

xf = x0;

xf.base_x = xstar.base_x + 1;
xf.base_y = xstar.base_y + 0.1;
xf.base_z = xstar.base_z;
nq = r.getNumPositions();

display('loading obstacles');

box_size = [0.4, 2.0, 0.05];

box_tops = [0, 0, 0.03, 0,3*pi/180;
            0.6, 0, 0.03, 0,-3*pi/180]';

safe_regions = iris.TerrainRegion.empty();

[A, b] = poly2lincon([0.8, 1.25, 0.8, 1.25], [-1, -1, 1, 1]);
[A, b] = convert_to_cspace(A, b);
safe_regions(1) = iris.TerrainRegion(A, b, [], [], [0;0;0], [0;0;1]);
[A, b] = poly2lincon([0.2, 0.4, 0.2, 0.4], [-1, -1, 1, 1]);
[A, b] = convert_to_cspace(A, b);
safe_regions(2) = iris.TerrainRegion(A, b, [], [], [0;0;0], [0;0;1]);

for j = 1:size(box_tops, 2)
  rpy = [box_tops(5,j);0;box_tops(4,j)];
  b = RigidBodyBox(box_size, box_tops(1:3,j) + [0;0;-box_size(3)/2], rpy);
  r = r.addGeometryToBody('world', b);
  offset = rpy2rotmat(rpy) * [0;0;box_size(3)/2];
  [A, b] = poly2lincon(box_tops(1,j) + offset(1) + [-0.22, 0.22, 0.22, -0.22],...
                       box_tops(2,j) + offset(2) + [-1, -1, 1, 1]);
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
lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_pose');
lcmgl.glColor3f(0,0,0);
f_height = xf.base_z - xstar.base_z;
lcmgl.glPushMatrix();
lcmgl.circle(xf.base_x, xf.base_y, f_height, 0.15);
lcmgl.glPopMatrix();
lcmgl.switchBuffers();

lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_rot');
lcmgl.glColor3f(0,0,0);
lcmgl.glTranslated(xf.base_x, xf.base_y, f_height);
axis = rpy2axis([0,0,xf.base_yaw]);
axis = axis([4,1:3]); 
axis(1) = 180 / pi * axis(1);
lcmgl.glRotated(axis(1), axis(2), axis(3), axis(4));
lcmgl.glPushMatrix();
len = 0.15;
lcmgl.glTranslated(len / 2, 0, 0);
lcmgl.drawArrow3d(len, 0.02, 0.02, 0.005);
lcmgl.glPopMatrix();
lcmgl.switchBuffers();

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
lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_pose');
lcmgl.glColor3f(0,0,0);
f_height = xf.base_z - xstar.base_z;
lcmgl.glPushMatrix();
lcmgl.circle(xf.base_x, xf.base_y, f_height, 0.15);
lcmgl.glPopMatrix();
lcmgl.switchBuffers();

lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'goal_rot');
lcmgl.glColor3f(0,0,0);
lcmgl.glTranslated(xf.base_x, xf.base_y, f_height);
axis = rpy2axis([0,0,xf.base_yaw]);
axis = axis([4,1:3]); 
axis(1) = 180 / pi * axis(1);
lcmgl.glRotated(axis(1), axis(2), axis(3), axis(4));
lcmgl.glPushMatrix();
len = 0.15;
lcmgl.glTranslated(len / 2, 0, 0);
lcmgl.drawArrow3d(len, 0.02, 0.02, 0.005);
lcmgl.glPopMatrix();
lcmgl.switchBuffers();