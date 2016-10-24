function [frame_knots, zmp_knots] = Quad_planSwingPitched(quadruped, stance1, stance2, stance3, swing1, swing2, initial_hold_time, target_frame_id)
% Compute a collision-free swing trajectory for a single foot.
if nargin < 5
  initial_hold_time = 0;
end

assert(swing1.frame_id == swing2.frame_id, 'planSwing expects to plan a swing trajectory between two positions of the /same/ foot body')
sizecheck(stance1.pos, [7, 1]);
sizecheck(stance2.pos, [7, 1]);
sizecheck(stance3.pos, [7, 1]);
sizecheck(swing1.pos, [7, 1]);
sizecheck(swing2.pos, [7, 1]);

params = struct(swing2.walking_params);
params = applyDefaults(params, quadruped.default_walking_params);
% NEW ----v
params.prevent_swing_overshoot = false;

DEBUG = false;
offset_feet = 0;
DEFAULT_FOOT_PITCH = 0; % Our quadruped has a single point end-effector

APEX_FRACTIONS = [0.05, 0.95]; % We plan only two poses of the foot during the aerial phase of the swing.
                               % Those poses are planned for locations where the toe has traveled a given
                               % fraction of the distance from its initial location to its final location.

FOOT_YAW_RATE = 0.375; % rad/s
MIN_STEP_TIME = 0.75; %s

MIN_DIST_FOR_PITCHED_SWING = 1.015;

if stance1.frame_id == quadruped.foot_frame_id.Foot1
  stance1_foot_name = 'Foot1';
elseif stance1.frame_id == quadruped.foot_frame_id.Foot2
  stance1_foot_name = 'Foot2';
elseif stance1.frame_id == quadruped.foot_frame_id.Foot3
  stance1_foot_name = 'Foot3';
else
  stance1_foot_name = 'Foot4';
end

if stance2.frame_id == quadruped.foot_frame_id.Foot1
  stance2_foot_name = 'Foot1';
elseif stance2.frame_id == quadruped.foot_frame_id.Foot2
  stance2_foot_name = 'Foot2';
elseif stance2.frame_id == quadruped.foot_frame_id.Foot3
  stance2_foot_name = 'Foot3';
else
  stance2_foot_name = 'Foot4';
end

if stance3.frame_id == quadruped.foot_frame_id.Foot1
  stance3_foot_name = 'Foot1';
elseif stance3.frame_id == quadruped.foot_frame_id.Foot2
  stance3_foot_name = 'Foot2';
elseif stance3.frame_id == quadruped.foot_frame_id.Foot3
  stance3_foot_name = 'Foot3';
else
  stance3_foot_name = 'Foot4';  
end

if swing1.frame_id == quadruped.foot_frame_id.Foot1
  swing_foot_name = 'Foot1';
elseif swing1.frame_id == quadruped.foot_frame_id.Foot2
  swing_foot_name = 'Foot2';
elseif swing1.frame_id == quadruped.foot_frame_id.Foot3
  swing_foot_name = 'Foot3';
else
  swing_foot_name = 'Foot4';
end

xy_dist = norm(swing2.pos(1:2) - swing1.pos(1:2));
% terrain_slice = double(swing2.terrain_pts);
% terrain_slice = [[0;swing1.pos(3)], terrain_slice, [xy_dist; swing2.pos(3)]];
% terrain_pts_in_local = [terrain_slice(1,:); zeros(1, size(terrain_slice, 2)); 
%                         terrain_slice(2,:) - swing1.pos(3)]

% Transform to world coordinates
T_local_to_world = [[rotmat(atan2(swing2.pos(2) - swing1.pos(2), swing2.pos(1) - swing1.pos(1))), [0;0];
                     0, 0, 1], [swing1.pos(1:2); 0]; 
                    0, 0, 0, 1];

% Determine how much of a forward step this is
R = quat2rotmat(stance1.pos(4:7));
swing_distance_in_local = (swing2.pos(1:3) - swing1.pos(1:3))' * (R * [1;0;0]);

toe_off_angle = 0;

swing_body_index = quadruped.foot_body_id.(swing_foot_name);
P1_body_index = quadruped.foot_body_id.(stance1_foot_name);
P2_body_index = quadruped.foot_body_id.(stance2_foot_name);
P3_body_index = quadruped.foot_body_id.(stance3_foot_name);
P4_body_index = quadruped.foot_body_id.(swing_foot_name);

% swing_toe_points_in_foot = quadruped.getBody(swing_body_index).getTerrainContactPoints('toe');
T_sole_to_foot = quadruped.getFrame(swing1.frame_id).T;
if target_frame_id.(swing_foot_name) < 0
  T_frame_to_foot = quadruped.getFrame(target_frame_id.(swing_foot_name)).T;
else
  T_frame_to_foot = eye(4);
end

T_swing1_sole_to_world = [quat2rotmat(swing1.pos(4:7)),swing1.pos(1:3); zeros(1, 3), 1];
T_swing1_frame_to_world = T_swing1_sole_to_world/T_sole_to_foot * T_frame_to_foot;

T_swing2_sole_to_world = [quat2rotmat(swing2.pos(4:7)),swing2.pos(1:3); zeros(1, 3), 1];
T_swing2_frame_to_world = T_swing2_sole_to_world/T_sole_to_foot * T_frame_to_foot;

toe1 = tform2poseQuat(T_swing1_frame_to_world/T_frame_to_foot);
toe2 = tform2poseQuat(T_swing2_frame_to_world/T_frame_to_foot);

quat_toe_off = rotmat2quat(T_swing1_frame_to_world(1:3,1:3) * rpy2rotmat([0;toe_off_angle;0]));
quat_swing2 = rotmat2quat(T_swing2_frame_to_world(1:3,1:3));

T_P1_sole_to_foot = quadruped.getFrame(stance1.frame_id).T;
if target_frame_id.(stance1_foot_name) < 0
  T_P1_frame_to_foot = quadruped.getFrame(target_frame_id.(stance1_foot_name)).T;
else
  T_P1_frame_to_foot = eye(4);
end
T_P1_sole_to_world = poseQuat2tform(stance1.pos);
T_P1_frame_to_world = T_P1_sole_to_world / T_P1_sole_to_foot * T_P1_frame_to_foot;

T_P2_sole_to_foot = quadruped.getFrame(stance2.frame_id).T;
if target_frame_id.(stance2_foot_name) < 0
  T_P2_frame_to_foot = quadruped.getFrame(target_frame_id.(stance2_foot_name)).T;
else
  T_P2_frame_to_foot = eye(4);
end
T_P2_sole_to_world = poseQuat2tform(stance2.pos);
T_P2_frame_to_world = T_P2_sole_to_world / T_P2_sole_to_foot * T_P2_frame_to_foot;

T_P3_sole_to_foot = quadruped.getFrame(stance3.frame_id).T;
if target_frame_id.(stance3_foot_name) < 0
  T_P3_frame_to_foot = quadruped.getFrame(target_frame_id.(stance3_foot_name)).T;
else
  T_P3_frame_to_foot = eye(4);
end
T_P3_sole_to_world = poseQuat2tform(stance3.pos);
T_P3_frame_to_world = T_P3_sole_to_world / T_P3_sole_to_foot * T_P3_frame_to_foot;

instep_shift = [0.0;stance1.walking_params.drake_instep_shift;0];

% NEW --v
zmp1 = (stance1.pos(1:2) + stance2.pos(1:2) + stance3.pos(1:2))/3;
%zmp1 = shift_step_inward(quadruped, stance1, instep_shift);

hold_time = params.drake_min_hold_time;

v1 = quat2rotmat(stance1.pos(4:7)) * [0;0;1];
b1 = -v1' * stance1.pos(1:3);
v2 = quat2rotmat(stance2.pos(4:7)) * [0;0;1];
b2 = -v2' * stance2.pos(1:3);
v3 = quat2rotmat(stance3.pos(4:7)) * [0;0;1];
b3 = -v3' * stance3.pos(1:3);

support_options.support_surface = {[v1;b1],[v2;b2],[v3;b3]};
support_options.contact_groups = {stance1.walking_params.support_contact_groups, stance2.walking_params.support_contact_groups, stance3.walking_params.support_contact_groups};
zmp_knots = struct('t', initial_hold_time + (hold_time / 2),...
 'zmp', zmp1, ...
 'supp', RigidBodySupportState(quadruped, [P1_body_index, P2_body_index, P3_body_index], support_options));

swing1_frame_pose = tform2poseQuat(T_swing1_frame_to_world);
swing2_frame_pose = tform2poseQuat(T_swing2_frame_to_world);

P1_frame_pose = tform2poseQuat(T_P1_frame_to_world);
P2_frame_pose = tform2poseQuat(T_P2_frame_to_world);
P3_frame_pose = tform2poseQuat(T_P3_frame_to_world);

frame_knots = struct('t', zmp_knots(end).t, ...
                    'Foot1', zeros(12,1),...
                    'Foot2', zeros(12,1),...
                    'Foot3', zeros(12,1),...
                    'Foot4', zeros(12,1),...
                    'toe_off_allowed', struct('Foot1', false, 'Foot2', false, 'Foot3', false, 'Foot4', false));

frame_knots.(stance1_foot_name) = [P1_frame_pose(1:3); quat2expmap(P1_frame_pose(4:7)); zeros(6,1)];
frame_knots.(stance2_foot_name) = [P2_frame_pose(1:3); quat2expmap(P2_frame_pose(4:7)); zeros(6,1)];
frame_knots.(stance3_foot_name) = [P3_frame_pose(1:3); quat2expmap(P3_frame_pose(4:7)); zeros(6,1)];

frame_knots.(swing_foot_name)  =  [swing1_frame_pose(1:3); quat2expmap(swing1_frame_pose(4:7)); zeros(6,1)];

function add_frame_knot(swing_pose, speed)
  if nargin < 2
    speed = params.step_speed/2;
  end
  frame_knots(end+1).(swing_foot_name) = [swing_pose(1:3); quat2expmap(swing_pose(4:7)); zeros(6,1)];
  frame_knots(end).(swing_foot_name)(4:6) = closestExpmap(frame_knots(end-1).(swing_foot_name)(4:6), frame_knots(end).(swing_foot_name)(4:6));
  frame_knots(end).(stance1_foot_name) = [P1_frame_pose(1:3); quat2expmap(P1_frame_pose(4:7)); zeros(6,1)];
  frame_knots(end).(stance3_foot_name) = [P3_frame_pose(1:3); quat2expmap(P3_frame_pose(4:7)); zeros(6,1)];
  frame_knots(end).(stance2_foot_name) = [P2_frame_pose(1:3); quat2expmap(P2_frame_pose(4:7)); zeros(6,1)];
  cartesian_distance = norm(frame_knots(end).(swing_foot_name)(1:3) - frame_knots(end-1).(swing_foot_name)(1:3));
  yaw_distance = abs((frame_knots(end).(swing_foot_name)(4:6) - frame_knots(end-1).(swing_foot_name)(4:6))' * [0;0;1]);
  dt = cartesian_distance / speed;
  frame_knots(end).t = frame_knots(end-1).t + dt;
  frame_knots(end).toe_off_allowed = struct('Foot1', false, 'Foot2', false, 'Foot3', false, 'Foot4', false);
end

max_terrain_ht_in_world = max(swing1.pos(3),swing2.pos(3));
% Apex knot 1
toe_apex1_in_world = (1-APEX_FRACTIONS(1))*toe1(1:3) + APEX_FRACTIONS(1)*toe2(1:3);
%initial_pose_not_flat = norm(cross(T_swing1_sole_to_world(1:3,1:3) * [0;0;1], [0;0;1])) >= sin(7.5*pi/180);
%if (initial_pose_not_flat || max_terrain_ht_in_world > toe_apex1_in_world(3) + params.step_height/4) && (~params.prevent_swing_undershoot)
%  toe_apex1_in_world = [toe1(1:2); max_terrain_ht_in_world + params.step_height];
%else
%  toe_apex1_in_world(3) = max([toe_apex1_in_world(3) + params.step_height,...
%                               max_terrain_ht_in_world + params.step_height]);
%end
toe_apex1_in_world(3) = max_terrain_ht_in_world + params.step_height;
T_apex1_toe_to_world = poseQuat2tform([toe_apex1_in_world(1:3); quat_toe_off]);
T_apex1_frame_to_world = T_apex1_toe_to_world * T_frame_to_foot;
T_apex1_frame_to_world(3,4) = toe_apex1_in_world(3);
add_frame_knot(tform2poseQuat(T_apex1_frame_to_world));

% Apex knot 2
toe_apex2_in_world = (1-APEX_FRACTIONS(2))*toe1(1:3) + APEX_FRACTIONS(2)*toe2(1:3);
%final_pose_not_flat = norm(cross(T_swing2_sole_to_world(1:3,1:3) * [0;0;1], [0;0;1])) >= sin(7.5*pi/180);
%if (final_pose_not_flat || max_terrain_ht_in_world > toe_apex2_in_world(3) + params.step_height/4) && (~params.prevent_swing_overshoot)
%  toe_apex2_in_world = [toe2(1:2); max_terrain_ht_in_world + params.step_height];
%else
%  toe_apex2_in_world(3) = max([toe_apex2_in_world(3) + params.step_height,...
%                               max_terrain_ht_in_world + params.step_height]);
%end
toe_apex2_in_world(3) = max_terrain_ht_in_world + params.step_height;
T_apex2_toe_to_world = poseQuat2tform([toe_apex2_in_world(1:3); quat_swing2]);
T_apex2_frame_to_world = T_apex2_toe_to_world * T_frame_to_foot;
T_apex2_frame_to_world(3,4) = toe_apex2_in_world(3);
add_frame_knot(tform2poseQuat(T_apex2_frame_to_world));

% Landing knot
add_frame_knot(swing2_frame_pose);
% add_foot_origin_knot(swing2_origin_pose, min(params.step_speed, MAX_LANDING_SPEED)/2);
if frame_knots(end).t - frame_knots(1).t < MIN_STEP_TIME
  frame_knots(end).t = frame_knots(1).t + MIN_STEP_TIME;
end

zmp_knots(end+1).t = frame_knots(end).t;
zmp_knots(end).zmp = zmp1;

vs1 = quat2rotmat(stance1.pos(4:7)) * [0;0;1];
bs1 = -v1' * stance1.pos(1:3);
vs2 = quat2rotmat(stance2.pos(4:7)) * [0;0;1];
bs2 = -v2' * stance2.pos(1:3);
vs3 = quat2rotmat(stance3.pos(4:7)) * [0;0;1];
bs3 = -v3' * stance3.pos(1:3);

v_swing = quat2rotmat(swing2.pos(4:7)) * [0;0;1];
b_swing = -v_swing' * swing2.pos(1:3);
support_options.support_surface = {[vs1;bs1], [vs2;bs2], [vs3;bs3], [v_swing;b_swing]};

support_options.contact_groups = {stance1.walking_params.support_contact_groups, stance2.walking_params.support_contact_groups, stance3.walking_params.support_contact_groups, swing2.walking_params.support_contact_groups};
zmp_knots(end).supp = RigidBodySupportState(quadruped, [P1_body_index, P2_body_index, P3_body_index, P4_body_index], support_options);

% Final knot
frame_knots(end+1) = frame_knots(end);
frame_knots(end).t = frame_knots(end-1).t + hold_time / 2;

% Find velocities for the apex knots by solving a small QP to get a smooth, minimum-acceleration cubic spline
foot = swing_foot_name;
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
  frame_knots(j).(foot)(7:12) = coefs(:,j,end-1);
end

end

function pos = shift_step_inward(quadruped, step, instep_shift)
  if step.frame_id == quadruped.foot_frame_id.Foot1 || step.frame_id == quadruped.foot_frame_id.Foot3
    instep_shift = [1;-1;1].*instep_shift;
  end
  if step.frame_id == quadruped.foot_frame_id.Foot1
    foot = 'Foot1';
  elseif step.frame_id == quadruped.foot_frame_id.Foot2
    foot = 'Foot2';
  elseif step.frame_id == quadruped.foot_frame_id.Foot3
    foot = 'Foot3';
  elseif step.frame_id == quadruped.foot_frame_id.Foot4
    foot = 'Foot4';
  else
    error('unknown frame ID: %d\n', step.frame_id);
  end
  T_sole_to_world = poseQuat2tform(step.pos);
  T_sole_to_foot = quadruped.getFrame(step.frame_id).T;

  T_foot_to_world = T_sole_to_world / T_sole_to_foot;

  supp_pts = quadruped.getTerrainContactPoints(quadruped.foot_body_id.(foot), {step.walking_params.support_contact_groups});
  supp_pts_in_world = T_foot_to_world * [supp_pts.pts; ones(1, size(supp_pts.pts, 2))];
  pose_center = [mean(supp_pts_in_world(1:3,:), 2); step.pos(4:7)];

  R = quat2rotmat(pose_center(4:7));
  shift = R*instep_shift;
  pos = pose_center(1:2) + shift(1:2);
end
