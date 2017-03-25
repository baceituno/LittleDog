function [zmp_knots, body_motions] = planSwing(quadruped, q0, footsteps, options)
% Plan the trajectories of the ZMP and the feet in order to follow the given footsteps
% @param q0 the initial configuration vector
% @param footsteps a list of Footstep objects
% @option t0 the initial time offset of the trajectories to be generated (default: 0)
% @option first_step_hold_s the number of seconds to wait before lifting the first foot (default: 1)

MAX_FINAL_SWING_SPEED = 1.0;

if nargin < 4; options = struct(); end

options = applyDefaults(options, struct('t0', 0,...
                                        'first_step_hold_s', 1.5,...
                                        'gait_sequence', repmat(eye(4),1,(length(footsteps)-4)/4)));

gait = options.gait_sequence;

target_frame_id = struct('Foot1', quadruped.foot_frame_id.Foot1,...
                         'Foot2', quadruped.foot_frame_id.Foot2,...
                         'Foot3', quadruped.foot_frame_id.Foot3,...
                         'Foot4', quadruped.foot_frame_id.Foot4);

typecheck(quadruped,{'RigidBodyManipulator','TimeSteppingRigidBodyManipulator'});
typecheck(q0,'numeric');
sizecheck(q0,[quadruped.getNumPositions,1]);

kinsol = doKinematics(quadruped, q0);
com0 = quadruped.getCOM(kinsol);
foot0 = struct('Foot1', forwardKin(quadruped, kinsol, quadruped.foot_frame_id.Foot1, [0;0;0], 2),...
               'Foot2', forwardKin(quadruped, kinsol, quadruped.foot_frame_id.Foot2, [0;0;0], 2),...
               'Foot3', forwardKin(quadruped, kinsol, quadruped.foot_frame_id.Foot3, [0;0;0], 2),...
               'Foot4', forwardKin(quadruped, kinsol, quadruped.foot_frame_id.Foot4, [0;0;0], 2));

footsteps_with_quat = footsteps;
footsteps_with_quat(1).pos = [footsteps_with_quat(1).pos(1:3); rpy2quat(footsteps_with_quat(1).pos(4:6))];

for j = 1:4
  for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
    foot = f{1};
    if footsteps_with_quat(j).frame_id == quadruped.foot_frame_id.(foot)
      footsteps_with_quat(j).pos = foot0.(foot);
    end
  end
end

for j = 1:length(footsteps)-1
  % Make sure quaternions of adjacent steps are in the same hemisphere
  N = 1;
  q1 = footsteps_with_quat(j).pos(4:7);
  q2 = rpy2quat(footsteps(j+N).pos(4:6));
  q2 = sign(q1' * q2) * q2;
  footsteps_with_quat(j+N).pos(4:7) = q2;
end

steps.Foot1 = footsteps_with_quat([footsteps_with_quat.frame_id] == quadruped.foot_frame_id.Foot1);
steps.Foot2 = footsteps_with_quat([footsteps_with_quat.frame_id] == quadruped.foot_frame_id.Foot2);
steps.Foot3 = footsteps_with_quat([footsteps_with_quat.frame_id] == quadruped.foot_frame_id.Foot3);
steps.Foot4 = footsteps_with_quat([footsteps_with_quat.frame_id] == quadruped.foot_frame_id.Foot4);

[steps.Foot1(1), steps.Foot2(1), steps.Foot3(1), steps.Foot4(1)] = getSafeInitialSupports(quadruped, kinsol, struct('Foot1', steps.Foot1(1), 'Foot2', steps.Foot2(1), 'Foot3', steps.Foot3(1), 'Foot4', steps.Foot4(1)));
[~, supp0] = getZMPBetweenFeet(quadruped, struct('Foot1', steps.Foot1(1), 'Foot2', steps.Foot2(1), 'Foot3', steps.Foot3(1), 'Foot4', steps.Foot4(1)));

% start zmp at current COM position
zmp0 = com0(1:2);

zmp_knots = struct('t', options.t0, 'zmp', zmp0, 'supp', supp0);

frame_knots = struct('t', options.t0, ...
  'Foot1', zeros(12,1),...
  'Foot2', zeros(12,1),...
  'Foot3', zeros(12,1),...
  'Foot4', zeros(12,1),...
  'toe_off_allowed', struct('Foot1', false, 'Foot2', false, 'Foot3', false, 'Foot4', false));

for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
  foot = f{1};
  frame_id = quadruped.foot_frame_id.(foot);
  T = quadruped.getFrame(frame_id).T;
  sole_pose = steps.(foot)(1).pos;
  Tsole = [quat2rotmat(sole_pose(4:7)), sole_pose(1:3); 0 0 0 1];
  Torig = Tsole / T;
  if target_frame_id.(foot) < 0
    Tframe = Torig * quadruped.getFrame(target_frame_id.(foot)).T;
  else
    Tframe = Torig;
  end
  frame_knots.(foot) = [Tframe(1:3,4); quat2expmap(rotmat2quat(Tframe(1:3,1:3))); zeros(6,1)];
end

istep = struct('Foot1', 1, 'Foot2', 1, 'Foot3', 1, 'Foot4', 1);
is_first_step = true;
cnt = 1;
while 1
  
  n_trans = gait(1,cnt) + gait(2,cnt) + gait(3,cnt) + gait(4,cnt);
  sw_feet = {};
  st_feet = {};

  n_sw = 0; n_st = 0;
  if(gait(1,cnt) == 1)
    sw_feet{n_sw + 1} = 'Foot1';
    n_sw = n_sw + 1;
  else
    st_feet{n_st + 1} = 'Foot1';
    n_st = n_st + 1;
  end

  if(gait(2,cnt) == 1)
    sw_feet{n_sw + 1} = 'Foot2';
    n_sw = n_sw + 1;
  else
    st_feet{n_st + 1} = 'Foot2';
    n_st = n_st + 1;
  end

  if(gait(3,cnt) == 1)
    sw_feet{n_sw + 1} = 'Foot3';
    n_sw = n_sw + 1;
  else
    st_feet{n_st + 1} = 'Foot3';
    n_st = n_st + 1;
  end

  if(gait(4,cnt) == 1)
    sw_feet{n_sw + 1} = 'Foot4';
    n_sw = n_sw + 1;
  else
    st_feet{n_st + 1} = 'Foot4';
    n_st = n_st + 1;
  end

  assert(n_st == 4-n_trans);
  assert(n_sw == n_trans);
  
  sw  = {};
  sw1 = {};
  st  = {};

  for l = 1:n_trans
    sw{l} = steps.(sw_feet{l})(istep.(sw_feet{l}));
    sw1{l} = steps.(sw_feet{l})(istep.(sw_feet{l})+1);
  end 

  for l = 1:(4-n_trans)
    st{l} = steps.(st_feet{l})(istep.(st_feet{l}));
  end

  if is_first_step
    initial_hold = options.first_step_hold_s;
    %sw1.walking_params.drake_min_hold_time = options.first_step_hold_s;
    if true
      is_first_step = false;
    % sw1.walking_params.step_speed = sw1.walking_params.step_speed / 2;
    end
  else
    initial_hold = 0;
  end

  if istep.Foot1 == length(steps.Foot1) || istep.Foot2 == length(steps.Foot2) || istep.Foot3 == length(steps.Foot3) || istep.Foot4 == length(steps.Foot4)
    % this is the last swing, so slow down
    for l = 1:n_trans
      sw1{l}.walking_params.step_speed = min(sw1{l}.walking_params.step_speed, MAX_FINAL_SWING_SPEED);
    end
  end

  [new_foot_knots, new_zmp_knots] = Quad_planSwingPitched(quadruped, st, sw, sw1, initial_hold, target_frame_id);
 
  t0 = frame_knots(end).t;
  for k = 1:length(new_foot_knots)
    new_foot_knots(k).t = new_foot_knots(k).t + t0;
  end
  for k = 1:length(new_zmp_knots)
    new_zmp_knots(k).t = new_zmp_knots(k).t + t0;
  end
  
  for l = 1:n_trans
    frame_knots(end).toe_off_allowed.(sw_feet{l}) = true;
  end

  frame_knots = [frame_knots, new_foot_knots];
  zmp_knots = [zmp_knots, new_zmp_knots];

  for l = 1:n_trans
    istep.(sw_feet{l}) = istep.(sw_feet{l}) + 1;
  end

  if istep.Foot1 == length(steps.Foot1) && istep.Foot2 == length(steps.Foot2) && istep.Foot3 == length(steps.Foot3) && istep.Foot4 == length(steps.Foot4)
    break
  end

  cnt = cnt + 1;
end

% add a segment at the end to recover
t0 = frame_knots(end).t;
frame_knots(end+1) = frame_knots(end);
frame_knots(end).t = t0 + 1.5;
[zmpf, suppf] = getZMPBetweenFeet(quadruped, struct('Foot1', steps.Foot1(end), 'Foot2', steps.Foot2(end), 'Foot3', steps.Foot3(end), 'Foot4', steps.Foot4(end)));
zmp_knots(end+1) =  struct('t', frame_knots(end-1).t, 'zmp', zmpf, 'supp', suppf);
zmp_knots(end+1) =  struct('t', frame_knots(end).t, 'zmp', zmpf, 'supp', suppf);

% dynamic_footstep_plan = DynamicFootstepPlan(quadruped, zmp_knots, frame_knots);
toe_off_allowed = [frame_knots.toe_off_allowed];
body_motions = BodyMotionData.empty();
for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
  foot = f{1};
  foot_states = [frame_knots.(foot)];
  frame_id = target_frame_id.(foot);
  ts = [frame_knots.t];
  body_motions(end+1) = BodyMotionData.from_body_xyzexp_and_xyzexpdot(frame_id, ts, foot_states(1:6,:), foot_states(7:12,:));
  body_motions(end).in_floating_base_nullspace = true(1, numel(ts));
  %% CHECK THE EFFECT OF THIS ---------------------v (maybe should be false)
  body_motions(end).toe_off_allowed = [toe_off_allowed.(foot)];
end

end

function [zmp, supp] = getZMPBetweenFeet(quadruped, steps)
  % Find the location of the ZMP at the center of the active contact points of the feet.
  % @param quadruped a quadruped
  % @param steps a struct mapping the strings to six Footstep objects
  zmp = zeros(2,0);
  initial_supports = [];
  initial_support_groups = {};
  initial_support_surfaces = {};
  for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
    foot = f{1};
    if steps.(foot).is_in_contact
      supp_groups = steps.(foot).walking_params.support_contact_groups;
      supp_pts = quadruped.getTerrainContactPoints(quadruped.foot_body_id.(foot), {supp_groups});
      initial_support_groups{end+1} = supp_groups;
      T_sole_to_world = poseQuat2tform(steps.(foot).pos);
      T_sole_to_foot = quadruped.getFrame(steps.(foot).frame_id).T;
      T_foot_to_world = T_sole_to_world / T_sole_to_foot;
      supp_pts_in_world = T_foot_to_world * [supp_pts.pts; ones(1, size(supp_pts.pts, 2))];
      zmp(:,end+1) = mean(supp_pts_in_world(1:2,:), 2);
      initial_supports(end+1) = quadruped.foot_body_id.(foot);
      v = quat2rotmat(steps.(foot).pos(4:7)) * [0;0;1];
      b = -v' * steps.(foot).pos(1:3);
      initial_support_surfaces{end+1} = [v;b];
    end
  end
  zmp = mean(zmp, 2);
  support_options.support_surface = initial_support_surfaces;
  support_options.contact_groups = initial_support_groups;
  supp = RigidBodySupportState(quadruped, initial_supports, support_options);
end

function [Foot1, Foot2, Foot3, Foot4] = getSafeInitialSupports(quadruped, kinsol, steps, options)
  % Make sure that our center of mass is within our initial available contact points by modifying the contact groups on the first two steps if necessary. 
  if nargin < 4
    options = struct('support_shrink_factor', 0.8);
  end
  com = getCOM(quadruped, kinsol);
  all_supp_pts_in_world = zeros(3, 0);
  for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
    foot = f{1};
    supp_groups = steps.(foot).walking_params.support_contact_groups;
    supp_pts = quadruped.getTerrainContactPoints(quadruped.foot_body_id.(foot), {supp_groups});
    supp_pts = [supp_pts.pts];
    supp_pts = bsxfun(@plus, mean(supp_pts, 2), options.support_shrink_factor * bsxfun(@minus, supp_pts, mean(supp_pts, 2)));
    supp_pts_in_world = quadruped.forwardKin(kinsol, quadruped.foot_body_id.(foot), supp_pts, 0);
    all_supp_pts_in_world = [all_supp_pts_in_world, supp_pts_in_world(1:3,:)];
  end
  k = convhull(all_supp_pts_in_world(1,:), all_supp_pts_in_world(2,:));
  if ~inpolygon(com(1), com(2), all_supp_pts_in_world(1,k), all_supp_pts_in_world(2,k))
    warning('Drake:CommandedSupportsDoNotIncludeCoM', 'Commanded support groups do not include the initial center of mass pose. Expanding the initial supports to prevent a fall at the start of walking');
    % CoM is outside the initial commanded support. This is almost certain to cause the robot to fall. So, for the first supports (the ones corresponding to the robot's current foot positions), we will allow the controller to use the entire heel-to-toe surface of the foot.
    for f = {'Foot1', 'Foot2', 'Foot3', 'Foot4'}
      foot = f{1};
      steps.(foot).walking_params.support_contact_groups = {'feet'};
    end
  end
  Foot1 = steps.Foot1;
  Foot2 = steps.Foot2;
  Foot3 = steps.Foot3;
  Foot4 = steps.Foot4;
end
