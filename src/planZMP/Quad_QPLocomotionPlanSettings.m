 classdef Quad_QPLocomotionPlanSettings
  properties
    robot;
    contact_groups;
    support_times
    supports;
    body_motions;
    zmp_data;
    zmptraj = [];
    V;
    qtraj;
    comtraj = [];
    mu = 0.7;
    D_control;
    use_plan_shift = false;
    plan_shift_body_motion_inds = 3;
    g = 9.81; % gravity m/s^2
    is_quasistatic = false;
    constrained_dofs = [];
    untracked_joint_inds;
    x0;

    %planned_support_command = Quad_QPControllerPlan.support_logic_maps.require_support; % when the plan says a given body is in support, require the controller to use that support. To allow the controller to use that support only if it thinks the body is in contact with the terrain, try QPControllerPlan.support_logic_maps.kinematic_or_sensed; 

    min_knee_angle = 0;
    knee_kp = 40;
    knee_kd = 4;
    knee_weight = 1;
    zmp_safety_margin = 0.005;

    % If a body is about to come into contact, then, after this fraction of the
    % duration of the current support, add an optional support that the
    % controller can use if it senses force on that body
    early_contact_allowed_fraction = 0.75;

    duration = inf;
    start_time = 0; 
    gain_set = 'standing';
  end

  methods
    function obj = Quad_QPLocomotionPlanSettings(robot)
      obj.robot = robot;
      %S = load(obj.robot.fixed_point_file);
      prop_cache = robot.getRobotPropertyCache();
      obj.contact_groups = prop_cache.contact_groups;
      q_fp = robot.getNomState();
      obj.qtraj = q_fp(1:prop_cache.nq);
      obj.constrained_dofs = [];
      obj.untracked_joint_inds = [];
      obj.zmp_data = Quad_QPLocomotionPlanSettings.defaultZMPData();
    end

    function obj = setCOMTraj(obj)
      ts = obj.qtraj.getBreaks();
      % this deals with the case of a constant trajectory
      com_poses = zeros(2,length(ts));
      for j = 1:numel(ts)
        kinsol = obj.robot.doKinematics(obj.qtraj.eval(ts(j)));
        com_position = obj.robot.getCOM(kinsol);
        com_poses(:,j) = com_position(1:2);
      end
      obj.comtraj = PPTrajectory(pchip(ts,com_poses));      
    end

    function obj = setLQRForCoM(obj)
      error('this has never been properly tuned on the robot and should not be used yet')
      Q = diag([10 10 1 1]);
      R = 0.0001*eye(2);
      A = [zeros(2),eye(2); zeros(2,4)];
      B = [zeros(2); eye(2)];
      [~,S,~] = lqr(A,B,Q,R);
      % set the Qy to zero since we only want to stabilize COM
      obj.zmp_data.Qy = 0*obj.zmp_data.Qy;
      obj.zmp_data.A = A;
      obj.zmp_data.B = B;
      obj.zmp_data.R = R;
      obj.zmp_data.S = S;
    end
    
   function draw_lcmgl(obj, lcmgl)    
      function plot_traj_foh(traj, color)    
        ts = traj.getBreaks();   
        pts = traj.eval(ts);  
        if size(pts,1) == 2    
          pts = [pts; zeros(1,size(pts,2))];   
        end    
        lcmgl.glColor3f(color(1), color(2), color(3));   
        lcmgl.glBegin(lcmgl.LCMGL_LINES);    
        for j = 1:length(ts)-1   
          lcmgl.glVertex3f(pts(1,j), pts(2,j),pts(3,j));   
          lcmgl.glVertex3f(pts(1,j+1), pts(2,j+1), pts(3,j+1));    
        end    
        lcmgl.glEnd();   
      end    
      link_trajectories = obj.getLinkTrajectories();   
      for j = 1:length(link_trajectories)    
        if ~isempty(link_trajectories(j).traj)
          plot_traj_foh(link_trajectories(j).traj, [0.8, 0.8, 0.2]);   
        else
          plot_traj_foh(link_trajectories(j).traj_min, [1, 0, 0]);   
          plot_traj_foh(link_trajectories(j).traj_max, [0, 1, 0]);   
        end    
      end    
      if ~isa(obj.comtraj, 'Trajectory')
        obj.comtraj = ExpPlusPPTrajectory(obj.comtraj.breaks,...   
                                          obj.comtraj.K,...    
                                          obj.comtraj.A,...    
                                          obj.comtraj.alpha,...    
                                          obj.comtraj.gamma);    
      end    
      %plot_traj_foh(obj.comtraj, [1,0,0]);   
      %plot_traj_foh(obj.zmptraj, [0,0,1]);   
    end    
   
    function link_trajectories = getLinkTrajectories(obj)    
      link_trajectories = struct('link_ndx', {}, 'traj', {}, 'min_traj', {}, 'max_traj', {});    
      for j = 1:length(obj.body_motions)   
        link_trajectories(j).link_ndx = obj.body_motions(j).body_id;   
        link_trajectories(j).traj = PPTrajectory(mkpp(obj.body_motions(j).ts, obj.body_motions(j).coefs, size(obj.body_motions(j).coefs, 1)));   
      end    
    end
    
  end

  methods(Static)
    
    function obj = fromQuadrupedFootstepPlan(footstep_plan, Quadruped, x0)
      zmp_options = struct();
      for j = 1:length(footstep_plan.footsteps)
        footstep_plan.footsteps(j).walking_params = applyDefaults(struct(footstep_plan.footsteps(j).walking_params),...
          Quadruped.default_walking_params);
      end
      
      [zmp_knots, foot_motion_data] = Quad_planZMPTraj(Quadruped,x0(1:Quadruped.getNumPositions()), footstep_plan.footsteps, zmp_options);
      obj = Quad_QPLocomotionPlanSettings.fromQuadrupedFootAndZMPKnots(foot_motion_data, zmp_knots, Quadruped, x0);
    end

    function obj = fromQuadrupedFootAndZMPKnots(foot_motion_data, zmp_knots, Quadruped, x0, options)
      if nargin < 5
        options = struct();
      end
      options = applyDefaults(options, struct('base_height_above_sole', Quadruped.default_walking_params.base_height_above_foot_sole,...
        'base_height_transition_knot',1));

      obj = Quad_QPLocomotionPlanSettings(Quadruped);
      obj.x0 = x0;
      % obj.qtraj = x0(1:Quadruped.getNumPositions());

      [obj.supports, obj.support_times] = Quad_QPLocomotionPlanSettings.getSupports(zmp_knots);
      obj.zmptraj = Quad_QPLocomotionPlanSettings.getZMPTraj(zmp_knots);
      [~, obj.V, obj.comtraj, LIP_height] = Quad_planZMPController(Quadruped,obj.zmptraj, obj.x0, options);
      obj.D_control = -Quadruped.default_walking_params.nominal_LIP_COM_height / obj.g * eye(2);
      obj.zmp_data.D = -LIP_height / obj.g * eye(2);
      obj.duration = obj.support_times(end)-obj.support_times(1)-0.001;
      if isa(obj.V.S, 'ConstantTrajectory')
        obj.V.S = fasteval(obj.V.S, 0);
      end
      obj.body_motions = [foot_motion_data];

      obj.gain_set = 'walking';
      obj.use_plan_shift = true;
    end

    function [supports, support_times] = getSupports(zmp_knots)
      supports = [zmp_knots.supp];
      support_times = [zmp_knots.t];
    end

    function zmptraj = getZMPTraj(zmp_knots)
      zmptraj = PPTrajectory(foh([zmp_knots.t], [zmp_knots.zmp]));
      zmptraj = setOutputFrame(zmptraj, SingletonCoordinateFrame('desiredZMP',2,'z',{'x_zmp','y_zmp'}));
    end
  end
  
  methods (Static, Access=private)
    function zmp_data = defaultZMPData()
      zmp_data = struct('A',  [zeros(2),eye(2); zeros(2,4)],... % COM state map 4x4
        'B', [zeros(2); eye(2)],... % COM input map 4x2
        'C', [eye(2),zeros(2)],... % ZMP state-output map 2x4
        'u0', zeros(2,1),... % nominal input 2x1
        'R', zeros(2),... % input LQR cost 2x2
        'Qy', eye(2));
    end
  end
end
