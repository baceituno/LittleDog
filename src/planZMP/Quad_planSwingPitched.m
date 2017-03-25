function [frame_knots, zmp_knots] = Quad_planSwingPitched(quadruped, st, sw, sw1, initial_hold_time, target_frame_id)
% Compute a collision-free swing trajectory for a timestep.
if nargin < 5
  initial_hold_time = 0;
end

nsw = size(sw,2);
nst = size(st,2);

for l = 1:nsw
  assert(sw{l}.frame_id == sw1{l}.frame_id, 'planSwing expects to plan a swing trajectory between two positions of the /same/ foot body')
  sizecheck(sw{l}.pos, [7, 1]);
  sizecheck(sw1{l}.pos, [7, 1]);
end

for l = 1:nst
  sizecheck(st{l}.pos, [7, 1]);
end

params = struct(sw{1}.walking_params);
params = applyDefaults(params, quadruped.default_walking_params);

DEBUG = false;
offset_feet = 0;
DEFAULT_FOOT_PITCH = 0; % Our quadruped has a single point end-effector

APEX_FRACTIONS = [0.05, 0.95]; % We plan only two poses of the foot during the aerial phase of the swing.
                               % Those poses are planned for locations where the toe has traveled a given
                               % fraction of the distance from its initial location to its final location.
MIN_STEP_TIME = 0.75; %s

stance_foot_name = {};
for l = 1:nst
  if st{l}.frame_id == quadruped.foot_frame_id.Foot1
    stance_foot_name{l} = 'Foot1';
  elseif st{l}.frame_id == quadruped.foot_frame_id.Foot2
    stance_foot_name{l} = 'Foot2';
  elseif st{l}.frame_id == quadruped.foot_frame_id.Foot3
    stance_foot_name{l} = 'Foot3';
  elseif st{l}.frame_id == quadruped.foot_frame_id.Foot4
    stance_foot_name{l} = 'Foot4';
  end
end

swing_foot_name = {};
for l = 1:nsw
  if sw{l}.frame_id == quadruped.foot_frame_id.Foot1
    swing_foot_name{l} = 'Foot1';
  elseif sw{l}.frame_id == quadruped.foot_frame_id.Foot2
    swing_foot_name{l} = 'Foot2';
  elseif sw{l}.frame_id == quadruped.foot_frame_id.Foot3
    swing_foot_name{l} = 'Foot3';
  elseif sw{l}.frame_id == quadruped.foot_frame_id.Foot4
    swing_foot_name{l} = 'Foot4';
  end
end

for l = 1:nsw
  xy_dist(l) = norm(sw1{l}.pos(1:2) - sw{l}.pos(1:2));

  % Transform to world coordinates
  T_local_to_world = [[rotmat(atan2(sw1{l}.pos(2) - sw{l}.pos(2), sw1{l}.pos(1) - sw{l}.pos(1))), [0;0];
                       0, 0, 1], [sw{l}.pos(1:2); 0]; 
                      0, 0, 0, 1];

  % Determine how much of a forward step this is
  R = quat2rotmat(st{l}.pos(4:7));
  swing_distance_in_local = (sw1{l}.pos(1:3) - sw{l}.pos(1:3))' * (R * [1;0;0]);
end

toe_off_angle = 0;
swing_body_index = {};
stance_body_index = {};

for l = 1:nsw
  swing_body_index{l} = quadruped.foot_body_id.(swing_foot_name{l});
end

for l = 1:nst
  stance_body_index{l} = quadruped.foot_body_id.(stance_foot_name{l});
end

T_swing1_sole_to_world = {};
T_swing1_frame_to_world = {};
T_swing2_sole_to_world = {};
T_swing2_frame_to_world = {};
T_frame_to_foot = {};
quat_toe_off = {};
quat_sw1 = {};
toe1 = {};
toe2 = {};
quat_toe_off = {};
quat_swing2 = {};

for l = 1:nsw

  T_sole_to_foot{l} = quadruped.getFrame(sw{l}.frame_id).T;
  if target_frame_id.(swing_foot_name{l}) < 0
    T_frame_to_foot{l} = quadruped.getFrame(target_frame_id.(swing_foot_name{l})).T;
  else
    T_frame_to_foot{l} = eye(4);
  end

  T_swing1_sole_to_world{l}  = [quat2rotmat(sw{l}.pos(4:7)),sw{l}.pos(1:3); zeros(1, 3), 1];
  T_swing1_frame_to_world{l} = T_swing1_sole_to_world{l}/T_sole_to_foot{l} * T_frame_to_foot{l};

  T_swing2_sole_to_world{l}  = [quat2rotmat(sw1{l}.pos(4:7)),sw1{l}.pos(1:3); zeros(1, 3), 1];
  T_swing2_frame_to_world{l} = T_swing2_sole_to_world{l}/T_sole_to_foot{l} * T_frame_to_foot{l};

  toe1{l} = tform2poseQuat(T_swing1_frame_to_world{l}/T_frame_to_foot{l});
  toe2{l} = tform2poseQuat(T_swing2_frame_to_world{l}/T_frame_to_foot{l});

  quat_toe_off{l} = rotmat2quat(T_swing1_frame_to_world{l}(1:3,1:3) * rpy2rotmat([0;toe_off_angle;0]));
  quat_swing2{l} = rotmat2quat(T_swing2_frame_to_world{l}(1:3,1:3));
end

T_P_sole_to_foot = {};
T_P_frame_to_foot = {};
T_P_sole_to_world = {};
T_P_frame_to_world = {};

for l = 1:nst
  T_P_sole_to_foot{l} = quadruped.getFrame(st{l}.frame_id).T;
  if target_frame_id.(stance_foot_name{l}) < 0
    T_P_frame_to_foot{l} = quadruped.getFrame(target_frame_id.(stance_foot_name{l})).T;
  else
    T_P_frame_to_foot{l} = eye(4);
  end
  T_P_sole_to_world{l} = poseQuat2tform(st{l}.pos);
  T_P_frame_to_world{l} = T_P_sole_to_world{l} / T_P_sole_to_foot{l} * T_P_frame_to_foot{l};
end

zmp1 = st{1}.pos(1:2);
if nst > 1
  for l = 2:nst
    zmp1 = zmp1 + st{l}.pos(1:2);
  end
end
zmp1(1:2) = zmp1(1:2)/nst;

hold_time = params.drake_min_hold_time;

v{l} = {};
b{l} = {};
for l = 1:nst
  v{l} = quat2rotmat(st{l}.pos(4:7)) * [0;0;1];
  b{l} = -v{l}' * st{l}.pos(1:3);
  support_options.support_surface{l} = [v{l};b{l}];
  support_options.contact_groups{l} = st{l}.walking_params.support_contact_groups;
end

arrayP_1 = [];
for l = 1:nst
  arrayP_1 = [arrayP_1, stance_body_index{l}];
end
zmp_knots = struct('t', initial_hold_time + (hold_time / 2),...
 'zmp', zmp1, ...
 'supp', RigidBodySupportState(quadruped, arrayP_1, support_options));

swing1_frame_pose = {};
swing2_frame_pose = {};
P_frame_pose = {};

for l = 1:nsw
  swing1_frame_pose{l} = tform2poseQuat(T_swing1_frame_to_world{l});
  swing2_frame_pose{l} = tform2poseQuat(T_swing2_frame_to_world{l});
end 

for l = 1:nst
  P_frame_pose{l} = tform2poseQuat(T_P_frame_to_world{l});
end 

frame_knots = struct('t', zmp_knots(end).t, ...
                    'Foot1', zeros(12,1),...
                    'Foot2', zeros(12,1),...
                    'Foot3', zeros(12,1),...
                    'Foot4', zeros(12,1),...
                    'toe_off_allowed', struct('Foot1', false, 'Foot2', false, 'Foot3', false, 'Foot4', false));

for l = 1:nsw
  frame_knots.(swing_foot_name{l})  =  [swing1_frame_pose{l}(1:3); quat2expmap(swing1_frame_pose{l}(4:7)); zeros(6,1)];
end 

for l = 1:nst
  frame_knots.(stance_foot_name{l}) = [P_frame_pose{l}(1:3); quat2expmap(P_frame_pose{l}(4:7)); zeros(6,1)];
end 

function add_frame_knot(swing_pose, speed)
  if nargin < 2
    speed = params.step_speed/2;
  end

  k = length(frame_knots);

  for l = 1:nsw
    frame_knots(k+1).(swing_foot_name{l}) = [swing_pose{l}(1:3); quat2expmap(swing_pose{l}(4:7)); zeros(6,1)];
    frame_knots(k+1).(swing_foot_name{l})(4:6) = closestExpmap(frame_knots(end-1).(swing_foot_name{l})(4:6), frame_knots(end).(swing_foot_name{l})(4:6));
  end

  for l = 1:nst
  frame_knots(end).(stance_foot_name{l}) = [P_frame_pose{l}(1:3); quat2expmap(P_frame_pose{l}(4:7)); zeros(6,1)];
  end

  cartesian_distance = norm(frame_knots(end).(swing_foot_name{1})(1:3) - frame_knots(end-1).(swing_foot_name{1})(1:3));
  yaw_distance = abs((frame_knots(end).(swing_foot_name{1})(4:6) - frame_knots(end-1).(swing_foot_name{1})(4:6))' * [0;0;1]);
  dt = cartesian_distance / speed;
  
  frame_knots(end).t = frame_knots(end-1).t + dt;
  frame_knots(end).toe_off_allowed = struct('Foot1', false, 'Foot2', false, 'Foot3', false, 'Foot4', false);
end

max_terrain_ht_in_world = cell(nsw);
toe_apex1_in_world = cell(nsw);
T_apex1_toe_to_world = cell(nsw);
T_apex1_frame_to_world = cell(nsw);
swing_pose = cell(nsw);

% Apex knot 1
for l = 1:nsw
  max_terrain_ht_in_world{l} = max(sw{l}.pos(3),sw1{l}.pos(3));
  toe_apex1_in_world{l} = (1-APEX_FRACTIONS(1))*toe1{l}(1:3) + APEX_FRACTIONS(1)*toe2{l}(1:3);
  toe_apex1_in_world{l}(3) = max_terrain_ht_in_world{l} + params.step_height;
  T_apex1_toe_to_world{l} = poseQuat2tform([toe_apex1_in_world{l}(1:3); quat_toe_off{l}]);
  T_apex1_frame_to_world{l} = T_apex1_toe_to_world{l} * T_frame_to_foot{l};
  T_apex1_frame_to_world{l}(3,4) = toe_apex1_in_world{l}(3);
  swing_pose{l} = tform2poseQuat(T_apex1_frame_to_world{l});
end
add_frame_knot(swing_pose);

toe_apex2_in_world = cell(nsw);
T_apex2_toe_to_world = cell(nsw);
T_apex2_frame_to_world = cell(nsw);
swing2_pose = cell(nsw);

% Apex knot 2
for l = 1:nsw
  toe_apex2_in_world{l} = (1-APEX_FRACTIONS(2))*toe1{l}(1:3) + APEX_FRACTIONS(2)*toe2{l}(1:3);
  toe_apex2_in_world{l}(3) = max_terrain_ht_in_world{l} + params.step_height;
  T_apex2_toe_to_world{l} = poseQuat2tform([toe_apex2_in_world{l}(1:3); quat_swing2{l}]);
  T_apex2_frame_to_world{l} = T_apex2_toe_to_world{l} * T_frame_to_foot{l};
  T_apex2_frame_to_world{l}(3,4) = toe_apex2_in_world{l}(3);
  swing2_pose{l} = tform2poseQuat(T_apex2_frame_to_world{l});
end
add_frame_knot(swing2_pose);
add_frame_knot(swing2_frame_pose);

if frame_knots(end).t - frame_knots(1).t < MIN_STEP_TIME
  frame_knots(end).t = frame_knots(1).t + MIN_STEP_TIME;
end

zmp_knots(end+1).t = frame_knots(end).t;
zmp_knots(end).zmp = zmp1;

vs = cell(nst);
bs = cell(nst);
v_swing = cell(nsw);
b_swing = cell(nsw);
for l = 1:nst
vs{l} = quat2rotmat(st{l}.pos(4:7)) * [0;0;1];
bs{l} = -v{l}' * st{l}.pos(1:3);
support_options.support_surface{l} = [vs{l};bs{l}];
support_options.contact_groups{l} = st{l}.walking_params.support_contact_groups;
end

for l = 1:nsw
  v_swing{l} = quat2rotmat(sw1{l}.pos(4:7)) * [0;0;1];
  b_swing{l} = -v_swing{l}' * sw1{l}.pos(1:3);
  support_options.support_surface{nst+l} = [v_swing{l};b_swing{l}];
  support_options.contact_groups{nst+l} = sw1{l}.walking_params.support_contact_groups;
end

arrayP_2 = [];

for l = 1:nst
  arrayP_2 = [arrayP_2, stance_body_index{l}];
end

for l = 1:nsw
  arrayP_2 = [arrayP_2, swing_body_index{l}];
end

zmp_knots(end).supp = RigidBodySupportState(quadruped, arrayP_2, support_options);

% Final knot
frame_knots(end+1) = frame_knots(end);
frame_knots(end).t = frame_knots(end-1).t + hold_time / 2;

% Find velocities for the apex knots by solving a small QP to get a smooth, minimum-acceleration cubic spline
foot = swing_foot_name{1};
states = [frame_knots(1:4).(foot)];
for j = 2:size(states, 2)
  % MAYBE USE N ----v
  w1 = states(4:6,j-1);
  w2 = states(4:6,j);
  [w2, dw2] = closestExpmap(w1, w2);
  states(4:6,j) = w2;
  states(10:12,j) = dw2 * states(10:12,j);
end

ts = [frame_knots(1).t, 0, 0, frame_knots(4).t];

qpSpline_options = struct('optimize_knot_times', true);
[coefs, ts] = qpSpline(ts,...
                 states(1:6,:),...
                 states(7:12, 1),...
                 states(7:12, 4), qpSpline_options);
%display('it passes :)');
for j = 2:3
  frame_knots(j).t = ts(j);
  for l = 1:nsw
    frame_knots(j).(swing_foot_name{1})(7:12) = coefs(:,j,end-1);
  end
end

end