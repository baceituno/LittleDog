classdef LittleDog < TimeSteppingRigidBodyManipulator & Quadruped
  methods

    function obj=LittleDog(urdf,options)
      typecheck(urdf,'char');

      if nargin < 2
        options = struct();
      end
      if ~isfield(options,'dt')
        options.dt = 0.001;
      end
      if ~isfield(options,'floating')
        options.floating = true;
      end
      if ~isfield(options,'terrain')
        options.terrain = RigidBodyFlatTerrain;
      end
      if ~isfield(options,'external_force')
        options.external_force = [];
      end
      
      w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');

      display('not yet')
      obj = obj@TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
      obj = obj@Quadruped('front_left_foot_center', 'front_right_foot_center', 'back_left_foot_center', 'back_right_foot_center');
      warning(w);

      % Add a force on a specified link if we want!
      if ~isempty(options.external_force)
        % For compile purposes, record that we have external force applied to a link
        % (this affects input frames)
        obj.external_force = findLinkId(obj,options.external_force);
        options_ef.weld_to_link = obj.external_force;
        obj = obj.addRobotFromURDF(fullfile(getDrakePath,'util','three_dof_force.urdf'), ...
          [0; 0; 0], [0; 0; 0], options_ef);
      end
      if options.floating
        % could also do fixed point search here
        obj = obj.setInitialState(obj.resolveConstraints(zeros(obj.getNumStates(),1)));
      else
        % TEMP HACK to get by resolveConstraints
        for i=1:length(obj.manip.body), obj.manip.body(i).contact_pts=[]; end
        obj.manip = compile(obj.manip);
        obj = obj.setInitialState(zeros(obj.getNumStates(),1));
      end
      
      obj.inner_foot_shape = [-0.01, -0.01, 0.01,  0.01;
                               0.01, -0.01, 0.01, -0.01];

      lastwarn = warning('off', 'Drake:RigidBodySupportState:NoSupportSurface');
      warning(lastwarn);
    end

    function obj = compile(obj)
      obj = compile@TimeSteppingRigidBodyManipulator(obj);
    end

    function obj = setInitialState(obj,x0)
      if isa(x0,'Point')
        obj.x0 = double(x0); %.inFrame(obj.getStateFrame));
      else
        typecheck(x0,'double');
        sizecheck(x0,obj.getNumStates());
        obj.x0 = x0;
      end
    end

    function x0 = getInitialState(obj)
      x0 = obj.x0;
    end

    function xstar = getNomState(obj)
      hip_roll = .1;
      hip_pitch = 1;
      knee = 1.55;

      xstar = Point(getStateFrame(obj));

      xstar.front_right_hip_roll = -hip_roll;
      xstar.front_right_hip_pitch = hip_pitch;
      xstar.front_right_knee = -knee;
      xstar.front_left_hip_roll = hip_roll;
      xstar.front_left_hip_pitch = hip_pitch;
      xstar.front_left_knee = -knee;
      xstar.back_right_hip_roll = -hip_roll;
      xstar.back_right_hip_pitch = -hip_pitch;
      xstar.back_right_knee = knee;
      xstar.back_left_hip_roll = hip_roll;
      xstar.back_left_hip_pitch = -hip_pitch;
      xstar.back_left_knee = knee;
      xstar.base_z = 0.1442;
      xstar.base_yaw = 0;
    end

    function weights = getFootstepOptimizationWeights(obj)
      % Return a reasonable set of default weights for the footstep planner
      % optimization. The weights describe the following quantities:
      % 'relative': the contribution to the cost function of the
      %             displacement from one step to the next
      % 'relative_final': the cost contribution of the displacement of the
      %                   displacement of the very last step (this can be
      %                   larger than the normal 'relative' cost in
      %                   order to encourage the feet to be close together
      %                   at the end of a plan)
      % 'goal': the cost contribution on the distances from the last two
      %         footsteps to their respective goal poses.
      % Each weight is a 6 element vector, describing the weights on
      % [x, y, z, roll, pitch, yaw]

      weights = struct('relative', [1;1;1;0;0;0.5],...
                       'relative_final', [10;10;10;0;0;2],...
                       'goal', [1000;1000;10;0;0;10]);
    end

    function prop_cache = getRobotPropertyCache(obj)
      % Functions like findLinkId, getTerrainContactPoints, etc. can be too slow to call
      % in the inner loop of our controller or planner, so we cache some useful information
      % at setup time. 
      prop_cache = struct('contact_groups', [],...
                          'body_ids', struct(),...
                          'position_indices', struct(),...
                          'actuated_indices', [],...
                          'nq', 0,...
                          'nv', 0,...
                          'num_bodies', 0);

      % getTerrainContactPoints is pretty expensive, so we'll just call it
      % for all the bodies and cache the results
      nbod = length(obj.getManipulator().body);
      contact_group_cache = cell(1, nbod);
      for j = 1:nbod
        contact_group_cache{j} = struct();
        for f = 1:length(obj.getBody(j).collision_geometry_group_names)
          name = obj.getBody(j).collision_geometry_group_names{f};
          if obj.getBody(j).robotnum == 1
            contact_group_cache{j}.(name) = obj.getBody(j).getTerrainContactPoints(name);
          end
        end
      end

      prop_cache.contact_groups = contact_group_cache;

      prop_cache.nq = obj.getNumPositions();
      prop_cache.nv = obj.getNumVelocities();
      prop_cache.num_bodies = length(obj.getManipulator().body);

      for b = {'body', 'front_left_lower_leg', 'front_right_lower_leg', 'back_left_lower_leg', 'back_right_lower_leg',...
               'front_left_upper_leg', 'front_right_upper_leg', 'back_left_upper_leg', 'back_right_upper_leg',...
               'front_left_hip', 'front_right_hip', 'back_left_hip', 'back_right_hip'}
        prop_cache.body_ids.(b{1}) = obj.findLinkId(b{1});
      end

      for j = {'front_left_hip_roll', 'front_left_hip_pitch', 'front_left_knee', 'front_right_hip_roll', 'front_right_hip_pitch', 'front_right_knee',...
               'back_left_hip_roll', 'back_left_hip_pitch', 'back_left_knee', 'back_right_hip_roll', 'back_right_hip_pitch', 'back_right_knee'}
        prop_cache.position_indices.(j{1}) = obj.findPositionIndices(j{1});
      end

      prop_cache.actuated_indices = obj.getActuatedJoints();
    end
    
    function footstep_plan = planFootstepsAdditive(obj, start_pos_or_q, goal_pos_or_q, safe_regions, iterations, N)

      footstep_plan = obj.planFootsteps(start_pos_or_q, goal_pos_or_q, safe_regions,struct('step_params', struct('max_num_steps', N + 4)));

      steps = struct();
      Foot1 = [footstep_plan.footsteps(end-3).pos(1), footstep_plan.footsteps(end-3).pos(2),... 
              footstep_plan.footsteps(end-3).pos(3), 0, 0, footstep_plan.footsteps(end-3).pos(6)]';
      Foot2 = [footstep_plan.footsteps(end-2).pos(1), footstep_plan.footsteps(end-2).pos(2),... 
              footstep_plan.footsteps(end-2).pos(3), 0, 0, footstep_plan.footsteps(end-2).pos(6)]';
      Foot3 = [footstep_plan.footsteps(end-1).pos(1), footstep_plan.footsteps(end-1).pos(2),... 
              footstep_plan.footsteps(end-1).pos(3), 0, 0, footstep_plan.footsteps(end-1).pos(6)]';
      Foot4 = [footstep_plan.footsteps(end).pos(1), footstep_plan.footsteps(end).pos(2),... 
              footstep_plan.footsteps(end).pos(3), 0, 0, footstep_plan.footsteps(end).pos(6)]';
      steps = struct('Foot1', Foot1, 'Foot2', Foot2, 'Foot3', Foot3, 'Foot4', Foot4);
      for i = 1:iterations
        clear footstep_plan1
        footstep_plan1 = obj.planFootsteps(steps, goal_pos_or_q, safe_regions,struct('step_params', struct('max_num_steps', N)));

        steps = struct();
        Foot1 = [footstep_plan1.footsteps(end-3).pos(1), footstep_plan1.footsteps(end-3).pos(2),... 
                footstep_plan1.footsteps(end-3).pos(3), 0, 0, footstep_plan1.footsteps(end-3).pos(6)]';
        Foot2 = [footstep_plan1.footsteps(end-2).pos(1), footstep_plan1.footsteps(end-2).pos(2),... 
                footstep_plan1.footsteps(end-2).pos(3), 0, 0, footstep_plan1.footsteps(end-2).pos(6)]';
        Foot3 = [footstep_plan1.footsteps(end-1).pos(1), footstep_plan1.footsteps(end-1).pos(2),... 
                footstep_plan1.footsteps(end-1).pos(3), 0, 0, footstep_plan1.footsteps(end-1).pos(6)]';
        Foot4 = [footstep_plan1.footsteps(end).pos(1), footstep_plan1.footsteps(end).pos(2),... 
                footstep_plan1.footsteps(end).pos(3), 0, 0, footstep_plan1.footsteps(end).pos(6)]';
        steps = struct('Foot1', Foot1, 'Foot2', Foot2, 'Foot3', Foot3, 'Foot4', Foot4);

        footstep_plan = footstep_plan.concat(footstep_plan1);
      end
    end

    function [plan, solvertime] = planFootsteps(obj, start_pos_or_q, goal_pos_or_q, safe_regions, options)
    % planFootsteps: find a set of reachable foot positions from the start to
    % the goal.
    % @param start_pos_or_q a struct with fields 'Footi' for each feet OR a configuration vector
    %                  start_pos.Footi is the 6 DOF initial pose
    %                  of the ith feet.
    % @param goal_pos a struct with fields 'Footi'.
    %                 goal_pos.Footi is the desired 6 DOF pose
    %                 of the ith foot, and likewise for each foot
    % @params safe_regions a list of planar polytopic regions into which the footstep
    %                      locations must be placed. Can be empty. If
    %                      safe_regions == [], then the footstep locations are
    %                      not constrained in this way. Otherwise, each
    %                      footstpe must fall within at least one of the
    %                      defined safe_regions. safe_regions should be a list
    %                      of objects or structures which have fields 'A',
    %                      'b', 'point', and 'normal' which define a region in
    %                      v = [x,y,z,yaw] as
    %                      A*v <= b AND normal'*v(1:3) == normal'*point
    % @param options a struct of options
    %
    % @option method_handle (default: @footstepPlanner.alternatingMIQP) the footstep planning
    %                method to use, expressed as a function handle
    % @option step_params (default: struct()) specific parameters for footstep
    %                     planning. Attributes set here overload those in
    %                     obj.default_footstep_params

      if nargin < 5; options = struct(); end
      if nargin < 4; safe_regions = []; end
      if ~isfield(options, 'method_handle'); options.method_handle = @ral2017; end
      if ~isfield(options, 'step_params'); options.step_params = struct(); end
      options.step_params = obj.applyDefaultFootstepParams(options.step_params);

      if isnumeric(start_pos_or_q)
        start_pos = obj.feetPosition(start_pos_or_q);
      else
        typecheck(start_pos_or_q, 'struct');
        start_pos = start_pos_or_q;
      end

      if isnumeric(goal_pos_or_q)
        goal_pos = obj.feetPosition(goal_pos_or_q);
      else
        typecheck(goal_pos_or_q, 'struct');
        goal_pos = goal_pos_or_q;
      end

      sizecheck(start_pos.Foot1, [6,1]);
      sizecheck(start_pos.Foot2, [6,1]);
      sizecheck(start_pos.Foot3, [6,1]);
      sizecheck(start_pos.Foot4, [6,1]);
      sizecheck(goal_pos.Foot1, [6,1]);
      sizecheck(goal_pos.Foot2, [6,1]);
      sizecheck(goal_pos.Foot3, [6,1]);
      sizecheck(goal_pos.Foot4, [6,1]);

      if isempty(safe_regions)
        terrain = obj.getTerrain();
        if isempty(terrain)
          % Use the plane of the feet to define the terrain
          n = rpy2rotmat(start_pos.Foot2(4:6)) * [0;0;1];
          pt = start_pos.Foot2(1:3);
        else
          [z, n] = terrain.getHeight(start_pos.Foot2(1:2));
          pt = [start_pos.Foot2(1:2); z];
        end
        safe_regions = [struct('A', zeros(0,3), 'b', zeros(0,1), 'point', pt, 'normal', n)];
      end
      for j = 1:length(safe_regions)
        sizecheck(safe_regions(j).A, [NaN, 3]);
        sizecheck(safe_regions(j).b, [size(safe_regions(j).A, 1), 1]);
        sizecheck(safe_regions(j).point, [3,1]);
        sizecheck(safe_regions(j).normal, [3,1]);
        for k = 1:size(safe_regions(j).A, 1)
          n = norm(safe_regions(j).A(k,:));
          safe_regions(j).A(k,:) = safe_regions(j).A(k,:) / n;
          safe_regions(j).b(k) = safe_regions(j).b(k) / n;
        end
      end
      % Always leads left
      foot1 = 'Foot1';
      foot2 = 'Foot2';
      foot3 = 'Foot3';
      foot4 = 'Foot4';
      
      plan = Quad_FootstepPlan.blank_plan(obj, options.step_params.max_num_steps,...
             [obj.foot_frame_id.(foot1), obj.foot_frame_id.(foot2), obj.foot_frame_id.(foot3),obj.foot_frame_id.(foot4)],...
             options.step_params, safe_regions);
      
      plan.footsteps(1).pos = start_pos.(foot1);
      plan.footsteps(2).pos = start_pos.(foot2);
      plan.footsteps(3).pos = start_pos.(foot3);
      plan.footsteps(4).pos = start_pos.(foot4);

      weights = getFootstepOptimizationWeights(obj);
      
      try

          [plan, solvertime] = options.method_handle(obj, plan, weights, goal_pos);
      catch e
        if strcmp(e.identifier, 'Drake:MixedIntegerConvexProgram:InfeasibleProblem')
          warning('The footstep planning problem is infeasible. Returning just the initial footstep poses');
          plan = plan.slice(1:6);
          solvertime = 0;
        else
          rethrow(e);
        end
      end
      
      end
  end
    
  properties (SetAccess = protected, GetAccess = public)
    x0
    external_force = 0; % if nonzero, body id where force is being exerted
  end

  properties
    fixed_point_file = fullfile('../model/BH3R_fp.mat');
    default_footstep_params = struct('nom_upward_step', 0.05,... % m
                                      'nom_downward_step', 0.05,...% m
                                      'max_num_steps', 67,...
                                      'min_num_steps', 7,...
                                      'leading_foot', 1); % i: Foot i
    default_walking_params = struct('step_speed', 0.5,... % speed of the swing foot (m/s)
                                    'step_height', 0.01,... % approximate clearance over terrain (m)
                                    'drake_min_hold_time', 0.0,... % minimum time in full support (s)
                                    'drake_instep_shift', 0.0,... % Distance to shift ZMP trajectory inward toward the instep from the center of the foot (m)
                                    'mu', 0.709,... % friction coefficient (Bullet standard)
                                    'constrain_full_foot_pose', true,... % whether to constrain the swing foot roll and pitch
                                    'base_height_above_foot_sole', 0.146,... % default base height when walking
                                    'support_contact_groups', {{'feet'}},... % which contact groups are available for support when walking
                                    'prevent_swing_undershoot', true,... % prevent the first phase of the swing from going backwards while moving to the first knot point
                                    'prevent_swing_overshoot', true,... % prevent the final phase of the swing from moving forward of the last knot point
                                    'nominal_LIP_COM_height', 0.146); % nominal height used to construct D_ls for our linear inverted pendulum model
    base_link = 'body';
    
    Foot1 = 'front_left_foot_center';
    Foot2 = 'front_right_foot_center';
    Foot3 = 'back_left_foot_center';
    Foot4 = 'back_right_foot_center';

    control_config_file = fullfile(fileparts(mfilename('fullpath')),'config.yaml');
  end
end
