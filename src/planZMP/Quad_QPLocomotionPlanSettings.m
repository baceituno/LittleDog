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
      options = struct();
      for j = 1:length(footstep_plan.footsteps)
        footstep_plan.footsteps(j).walking_params = applyDefaults(struct(footstep_plan.footsteps(j).walking_params),...
          Quadruped.default_walking_params);
      end
      options.gait_sequence = footstep_plan.gait;
      [zmp_knots, foot_motion_data] = Quad_planSwing(Quadruped,x0(1:Quadruped.getNumPositions()), footstep_plan.footsteps, options);
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
