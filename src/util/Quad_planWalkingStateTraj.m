function [xtraj, htraj, ts] = Quad_planWalkingStateTraj(robot, walking_plan_data, xstar,comtraj)
% Given the results of the ZMP tracker, find a state trajectory for the robot to execute
% its walking plan.
% @param walking_plan_data a WalkingPlanData, such as that returned by Quad_planWalkingZMP()
% @param xstar the nominal robot state vector
% @retval xtraj a PPTrajectory of robot states
% @retval htraj a PPTrajectory of CoM heights
% @retval ts the time points at which the trajectory constraints were applied

xstar = robot.getNomState();

nq = robot.getNumPositions();
q0 = walking_plan_data.x0(1:nq);
qstar = xstar(1:nq);

% time spacing of samples for IK
if ~isa(walking_plan_data.comtraj, 'Trajectory')
  walking_plan_data.comtraj = ExpPlusPPTrajectory(walking_plan_data.comtraj.breaks, walking_plan_data.comtraj.K, walking_plan_data.comtraj.A, walking_plan_data.comtraj.alpha, walking_plan_data.comtraj.gamma);
end

ts = 0:0.1:walking_plan_data.comtraj.tspan(end);
if length(ts)>150 % limit number of IK samples to something reasonable
  ts = linspace(0,walking_plan_data.comtraj.tspan(end),150);
end

cws = 1;
if nargin < 4
  comtraj = walking_plan_data.comtraj;
  cws = 0;
end

% create desired joint trajectory
cost = Point(robot.getStateFrame,1);
cost.base_roll = 1000;
cost.base_pitch = 1000;
cost.base_yaw = 0;
cost = 1.0 + double(cost);

q = zeros(robot.getNumPositions(), length(ts));
htraj = [];

length(ts)
for i = 1:length(ts)
  t = ts(i);
  i
  if (i > 1)
    ik_args = {};
    for j = 1:length(walking_plan_data.body_motions)
      body_ind = walking_plan_data.body_motions(j).body_id;
      xyz_exp = walking_plan_data.body_motions(j).eval(t);
      xyz = xyz_exp(1:3);
      quat = expmap2quat(xyz_exp(4:6));
      xyz(walking_plan_data.body_motions(j).weight_multiplier(4:6) == 0) = nan;
      ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,false,robot.getManipulator(),body_ind, [0;0;0],xyz, xyz)}];
    end
    if cws
      com = comtraj.eval(t);
      ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldCoMConstraintType ,false,robot.getManipulator(),[com(1)-0.04;com(2)-0.04;com(3)-0.02],[com(1)+0.04;com(2)+0.04;com(3)+0.02])}];
    end
    ik_prob = InverseKinematics(robot, qstar, ik_args{:});
    ik_prob = ik_prob.setQ(diag(cost(1:robot.getNumPositions)));
    ik_prob = ik_prob.setSolver('snopt');
    q(:,i) = ik_prob.solve(q(:,i-1));
    %q(:,i) = inverseKin(robot,q(:,i-1),qstar,Quad_joint_cnstr,kc_com,ik_args{:},ikoptions);
  else
    q = q0;
  end
  com = getCOM(robot,q(:,i));
  htraj = [htraj com(3)];
end

% qtraj = PPTrajectory(spline(ts,q));
htraj = PPTrajectory(spline(ts,htraj));
x = zeros(getNumStates(robot),length(ts));
x(1:getNumPositions(robot),:) = q;
xtraj = PPTrajectory(spline(ts, x));
xtraj = xtraj.setOutputFrame(robot.getStateFrame());

end


