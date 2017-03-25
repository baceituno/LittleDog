%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernardo Aceituno C.         %
% USB C Laboratory             %
% Mechatronics Research Group  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run LittleDog across a rough terrain using the Robust Walking Motion Planner

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

xf = x0;

xf.base_x = xstar.base_x + 3.0;
nq = r.getNumPositions();

safe_regions = iris.TerrainRegion.empty();

[A, b] = poly2lincon([-2, 8, 8, -2], [-1, -1, 1, 1]);
[A, b] = convert_to_cspace(A, b);
safe_regions(1) = iris.TerrainRegion(A, b, [], [], [0;0;0], [0;0;1]);
r = compile(r);

display('setting terrain');
height_map = RigidBodyHeightMapTerrain.constructHeightMapFromRaycast(r,x0(1:nq),-1:.015:4,-1:.015:1,10);
r = r.setTerrain(height_map).compile();
r.default_walking_params.drake_min_hold_time = 1.0;

% display('Building visualizer');
% v = r.constructVisualizer();
% v.display_dt = 0.01;
% v.draw(0,x0(1:nq))

display('planning footsteps')
footstep_plan = r.planFootstepsAdditive(x0(1:nq), xf(1:nq), safe_regions, 0, 12);

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
% lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'footstep_plan');
%   footstep_plan.draw_lcmgl(lcmgl);
%   lcmgl.switchBuffers();
display('Planning Body Motions')
settings = Quad_QPLocomotionPlanSettings.fromQuadrupedFootstepPlan(footstep_plan, r, x0(1:nq));

% display('Drawing plan')
% if isa(v, 'BotVisualizer')
%   lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'footstep_plan');
%   footstep_plan.draw_lcmgl(lcmgl);
%   lcmgl.switchBuffers();
%   lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'walking_plan');
%   settings.draw_lcmgl(lcmgl);
%   lcmgl.switchBuffers();
% else
%   figure(25)
%   footstep_plan.draw_2d();
% end

dt = 0.1;
mu_ground = 1.0;
num_fc_edges = 4;
s = 1/0.1331;
Q_w = diag([1 1 1 s s s]);
com = getCOM(r,xstar(1:nq));
com_height_guess = com(3);

centroidal_angular_momentum_norm_weight = [10;10;10];

display('parsing walking plan');
walking_plan_parser = Quad_WalkingPlanParser(footstep_plan,settings,dt,mu_ground,num_fc_edges,com_height_guess,Q_w);
robot_mass = r.getMass();

% ellipsoid com region
A_com = diag([0.01;0.01;0.005]);

display('declaring the CWS planner');
cws_planner = iros2017(r.getMass(), walking_plan_parser, Q_w, A_com, centroidal_angular_momentum_norm_weight);
cws_planner.vars.cws_margin.lb = 0*ones(walking_plan_parser.nT,1);
% initial state constraint
cws_planner.vars.r.lb(1:2,1) = walking_plan_parser.com_guess(1:2,1);
cws_planner.vars.r.ub(1:2,1) = walking_plan_parser.com_guess(1:2,1);
% final state constraint
cws_planner.vars.r.lb(1:2,end) = walking_plan_parser.com_guess(1:2,end);
cws_planner.vars.r.ub(1:2,end) = walking_plan_parser.com_guess(1:2,end);

display('solving for COM trajectory');
tic;[cws_planner,~,~] = cws_planner.solveMosek();toc;

com_traj = cws_planner.vars.r.value;

% lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(),'com_traj');
% lcmgl.glLineWidth(3);
% lcmgl.glColor3f(0,0,1);
% lcmgl.plot3(com_traj(1,:),com_traj(2,:),com_traj(3,:));
% lcmgl.switchBuffers();

%{
figure(1)
t = 1:length(kG_norm(idx_init:end));
plot(t, kG_norm(idx_init:end), 'b', t, kG_ub_norm(idx_init:end),  'r--')
xlim([0 t(end)])
ylabel({'$|k_G|_1$'},'interpreter','LaTeX','FontSize',12);
xlabel('$timestep$','interpreter','LaTeX','FontSize',12);
legend('|k_G|_1','s','Location','northeast')
print -depsc2 flat_kg.eps
%}

display('saving_vars ')
idx_init = 130;
kG = cws_planner.vars.kO.value - cws_planner.robot_mass*cross(cws_planner.vars.r.value,cws_planner.vars.dr.value,1);
kG = kG';
kG_norm = abs(kG(:,1)) + abs(kG(:,2)) + abs(kG(:,3));
kG_ub_norm = cws_planner.vars.kG_norm_ub.value(:);
kG_ub_norm = kG_ub_norm';

cws_m = cws_planner.vars.cws_margin.value();
cws_m = cws_m';
com_guess = walking_plan_parser.com_guess(:,idx_init:end);
com_plan = com_traj(:,idx_init:end);
figure(1)
subplot(2,1,1)
plot(com_traj(1,:))
subplot(2,1,2)
plot(com_traj(2,:))
figure(2)
plot(cws_m)

%running IK and simulating
ts = 0:0.1:settings.comtraj.tspan(end);
int = 1 + size(com_traj,2) - length(ts);
comtraj = PPTrajectory(spline(ts, com_traj(:,int:end)));

%[xtraj, htraj, ts] = Quad_planWalkingStateTraj(r,settings,xstar,comtraj);
%v.playback(xtraj, struct('slider', true));