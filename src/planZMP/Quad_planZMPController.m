function [c, V, comtraj, limp_height] = Hex_planZMPController(obj, zmptraj_or_comgoal, x0, options)
if nargin < 4
  options = struct();
end
options = applyDefaults(options, struct('base_height_above_sole', []));

path_handle = addpathTemporary(fullfile(getDrakePath(), 'examples', 'ZMP'));

nq = obj.getNumPositions();
q0 = x0(1:nq);

kinsol = doKinematics(obj, q0);

zmp_opts = struct('use_lqr_cache', true, 'lqr_cache_com_height_resolution', 0.01, 'build_control_objects', false);
if size(x0, 1) == nq
  com = getCOM(obj, kinsol);
  zmp_opts.com0 = com(1:2);
else
  assert(size(x0, 1) == obj.getNumStates(), 'expected a full configuration or state vector');
  [com, J] = getCOM(obj, kinsol);
  comdot = J * x0((obj.getNumPositions()+1):end);
  zmp_opts.com0 = com(1:2);
  zmp_opts.comdot0 = comdot(1:2);
end

foot_pos = [forwardKin(obj, kinsol, obj.foot_frame_id.Foot1, zeros(3, 1)), forwardKin(obj, kinsol, obj.foot_frame_id.Foot2, zeros(3, 1)),...
            forwardKin(obj, kinsol, obj.foot_frame_id.Foot3, zeros(3, 1)), forwardKin(obj, kinsol, obj.foot_frame_id.Foot4, zeros(3, 1))];

zfeet = mean(foot_pos(3,:));

base_pos = forwardKin(obj, kinsol, obj.findLinkId('body'), [0;0;0]);
actual_base_above_sole = base_pos(3) - zfeet;

desired_base_above_sole = options.base_height_above_sole;

limp_height = (com(3) - zfeet) + (desired_base_above_sole - actual_base_above_sole);
limp_height = round(limp_height / zmp_opts.lqr_cache_com_height_resolution) * zmp_opts.lqr_cache_com_height_resolution;
if isa(zmptraj_or_comgoal, 'Trajectory')
  [c, V, comtraj] = LinearInvertedPendulum.ZMPtrackerClosedForm(limp_height, zmptraj_or_comgoal, zmp_opts);
elseif isnumeric(zmptraj_or_comgoal)
  hg = limp_height/9.81;
  [~, S] = ZMPCachedLQR(hg, diag([0,0,0,0,1,1]));
  c = [];
  V = struct('S', S, 's1', zeros(4,1));
  comtraj = zmptraj_or_comgoal;
else
  error('Wrong type. Expected a zmp trajectory or a desired com position in xy');
end
