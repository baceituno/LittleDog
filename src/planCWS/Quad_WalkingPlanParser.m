classdef Quad_WalkingPlanParser
  % parse the footstep_plan and walking plan, so that we can construct a
  % ContactWrenchSetPlanner
  properties(SetAccess = protected)
    num_footsteps
    footstep_contact_pos
    support_times
    footstep_support_interval % A 2 x num_footsteps matrix, the i'th footstep starts in support_times(footstep_support_interval(1,i)) and ends at support_times(footstep_support_interval(2,i))
    
    footstep_plan
    walking_plan
    
    t
    nT
    robot
    robot_mass
    gravity
    wrench_disturbance_pt
    Ain_cws
    bin_cws
    com_guess;
    Linearized_cones;
    
    Q_w
    
    mu_ground
    num_fc_edges
    footstep_fc_edges
    num_supports
  end
  
  methods
    function obj = Quad_WalkingPlanParser(footstep_plan,walking_plan,dt,mu_ground,num_fc_edges,com_height_guess,Q_w)
      if(~isa(footstep_plan,'Quad_FootstepPlan') || ~isa(walking_plan,'Quad_QPLocomotionPlanSettings'))
        error('input should be a Quad_FootstepPlan and a Quad_QPLocomotionPlanSettings');
      end
      obj.footstep_plan = footstep_plan;
      obj.walking_plan = walking_plan;
      
      obj.robot = walking_plan.robot;
      obj.robot_mass = obj.robot.getMass();
      obj.gravity = norm(obj.robot.getGravity);
      obj.support_times = walking_plan.support_times;
      obj = obj.getSampleTime(dt);
      obj.mu_ground = mu_ground;
      obj.num_fc_edges = num_fc_edges;
      obj.Q_w = Q_w;
      obj = obj.getFootsteps();
      
      obj = obj.computeContactWrenchSet(com_height_guess);
      
      obj.wrench_disturbance_pt(:,end) = obj.wrench_disturbance_pt(:,end-1);
      obj.com_guess(:,end) = obj.com_guess(:,end-1);
      obj.Ain_cws{obj.nT} = obj.Ain_cws{obj.nT-1};
      obj.bin_cws{obj.nT} = obj.bin_cws{obj.nT-1};
    end
  end
  
  methods(Access = protected)
    function obj = getSampleTime(obj,dt)
      if(numel(dt)~=1 || dt<=0)
        error('dt should be a positive scalar');
      end
      start_time = obj.walking_plan.start_time;
      obj.t = start_time;
      for i = 1:length(obj.walking_plan.supports)-1
        t_new = obj.t(end):dt:obj.walking_plan.support_times(i+1);
        if(t_new(end) ~= obj.walking_plan.support_times(i+1))
          t_new = [t_new obj.walking_plan.support_times(i+1)];
        end
        obj.t = [obj.t;t_new(2:end)'];       
      end
      obj.nT = length(obj.t);
    end
    
    function edges = frictionConeEdges0(obj,mu_fc,num_edges)
      theta = linspace(0,2*pi,num_edges+1);
      theta = theta(1:end-1);
      edges = [mu_fc*cos(theta);mu_fc*sin(theta);ones(1,num_edges)];
    end
    
    function obj = getFootsteps(obj)
      fc_edges0 = obj.frictionConeEdges0(obj.mu_ground,obj.num_fc_edges);
      % Initially all feet on ground
      foot_ids = obj.walking_plan.supports(1).bodies; 
      if(length(foot_ids) ~= 4)
        error('Why there are not four feet on the ground?');
      end
      
      obj.num_footsteps = length(obj.footstep_plan.footsteps);
      % footstep_plan.footsteps(i) is active from
      % support_times(obj.footstep_support_interval(1,i)) to
      % support_times(obj.footstep_support_interval(2,i))
      obj.footstep_support_interval = zeros(2,obj.num_footsteps);
      obj.num_supports = length(obj.walking_plan.supports);
      obj.footstep_contact_pos = cell(obj.num_footsteps,1);
      obj.footstep_fc_edges = cell(obj.num_footsteps,1);
      for i = 1:length(obj.footstep_plan.footsteps) 
        T_sole_to_world = [rpy2rotmat(obj.footstep_plan.footsteps(i).pos(4:6)) obj.footstep_plan.footsteps(i).pos(1:3);zeros(1,3) 1];
        frame_id = obj.footstep_plan.footsteps(i).frame_id;
        T_foot_to_world = T_sole_to_world/obj.robot.getFrame(frame_id).T;
        body_id = obj.robot.getFrame(frame_id).body_ind;
        % First find the starting support time
        for j = obj.footstep_support_interval(1,i)+1:obj.num_supports
          if(any(obj.walking_plan.supports(j).bodies == body_id))
            obj.footstep_support_interval(1,i) = j;
            % Now find the ending support time
            find_end_time = false;
            for k = j+1:obj.num_supports
              if(~any(obj.walking_plan.supports(k).bodies == body_id))
                obj.footstep_support_interval(2,i) = k;
                find_end_time = true;
                break;
              end
            end
            if(~find_end_time)
              obj.footstep_support_interval(2,i) = obj.num_supports;
            end
            % For all subsequent footsteps with the same foot, they cannot
            % start earlier than the end time of the current foot
            for k = i+1:obj.num_footsteps
              if(obj.footstep_plan.footsteps(k).frame_id == frame_id)
                obj.footstep_support_interval(1,k) = obj.footstep_support_interval(2,i);
              end
            end
            break;
          end
        end
        % compute the contact point position
        body_in_support_idx = body_id == obj.walking_plan.supports(obj.footstep_support_interval(1,i)).bodies;
        body_pts = obj.walking_plan.supports(obj.footstep_support_interval(1,i)).contact_pts{body_in_support_idx};
        contact_pt_pos = [eye(3) zeros(3,1)]*T_foot_to_world*[body_pts;ones(1,size(body_pts,2))];
        obj.footstep_contact_pos{i} = contact_pt_pos;
        footstep_normal = obj.walking_plan.supports(obj.footstep_support_interval(1,i)).support_surface{body_in_support_idx}(1:3);
        R_fc = rotateVectorToAlign([0;0;1],footstep_normal);
        obj.footstep_fc_edges{i} = R_fc*fc_edges0;
        obj.Linearized_cones{i} = LinearizedFrictionCone(obj.footstep_contact_pos{i},footstep_normal,obj.mu_ground,obj.footstep_fc_edges{i});
      end
    end
    
    function obj = computeContactWrenchSet(obj,com_height_guess)
      % we apply the disturbance at the center of the foot region
      obj.wrench_disturbance_pt = zeros(3,obj.nT);
      % com_guess is a guess of the com trajectory, we then add a box
      % around the com_guess to bound the com position
      obj.com_guess = zeros(3,obj.nT);
      % Possible case
      % 1. contact feet do not change
      % 2. From not in contact to contact
      % 3. From contact to not in contact
      obj.Ain_cws = cell(obj.nT,1);
      obj.bin_cws = cell(obj.nT,1);
      %obj.cws_central = zeros(6,obj.nT);
      %obj.cws_central_margin = zeros(1,obj.nT);
      obj.num_supports
      for i = 2:obj.num_supports
        bodies0 = obj.walking_plan.supports(i-1).bodies;
        bodies1 = obj.walking_plan.supports(i).bodies;
        t0 = obj.walking_plan.support_times(i-1);
        t1 = obj.walking_plan.support_times(i);
        t_idx = find(obj.t>=t0 &obj.t<t1);
        if(length(bodies0) == 3 && length(bodies1)== 4)
          % from not in contact to contact
          % find the swing foot pos at the begining of the step
          swing_pos0 = obj.footstep_contact_pos{obj.footstep_support_interval(2,:) == i-1};
          swing_pos1 = obj.footstep_contact_pos{obj.footstep_support_interval(1,:) == i};
          stance_pos = obj.footstep_contact_pos(obj.footstep_support_interval(1,:)<=i-1&obj.footstep_support_interval(2,:)>i);
          % disturbance is at the center of the support surface
          obj.wrench_disturbance_pt(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2)],2));
          com_guess0 = mean([mean(swing_pos0,2) mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2)],2);
          com_guess0(3) = com_guess0(3)+com_height_guess;
          com_guess1 = mean([mean(swing_pos1,2) mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2)],2);
          com_guess1(3) = com_guess1(3)+com_height_guess;
          obj.com_guess(:,t_idx) = (com_guess1-com_guess0)/(t1-t0)*(obj.t(t_idx)-t0)'+repmat(com_guess0,1,length(t_idx));
        elseif(length(bodies0) == 2 && length(bodies1)== 4)
          % from not in contact to contact
          % find the swing foot pos at the begining of the step
          swing_pos0 = obj.footstep_contact_pos(obj.footstep_support_interval(2,:) == i-1);
          swing_pos1 = obj.footstep_contact_pos(obj.footstep_support_interval(1,:) == i);
          stance_pos = obj.footstep_contact_pos(obj.footstep_support_interval(1,:)<=i-1&obj.footstep_support_interval(2,:)>i);
          % disturbance is at the center of the support surface
          obj.wrench_disturbance_pt(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2)],2));
          com_guess0 = mean([mean(swing_pos0{1},2) mean(swing_pos0{2},2) mean(stance_pos{1},2) mean(stance_pos{2},2)],2);
          com_guess0(3) = com_guess0(3)+com_height_guess;
          com_guess1 = mean([mean(swing_pos1{1},2) mean(swing_pos1{2},2) mean(stance_pos{1},2) mean(stance_pos{2},2)],2);
          com_guess1(3) = com_guess1(3)+com_height_guess;
          obj.com_guess(:,t_idx) = (com_guess1-com_guess0)/(t1-t0)*(obj.t(t_idx)-t0)'+repmat(com_guess0,1,length(t_idx));
        elseif(length(bodies0) == 4 && length(bodies1) == 4)
          % both feet on the ground
          stance_pos = obj.footstep_contact_pos(obj.footstep_support_interval(1,:)<=i-1&obj.footstep_support_interval(2,:)>=i);
          % disturbance is at the center of the supports
          obj.wrench_disturbance_pt(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2) mean(stance_pos{4},2)],2));
          obj.com_guess(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2) mean(stance_pos{4},2)],2)+[0;0;com_height_guess]);
        elseif(length(bodies0) == 4 && length(bodies1) == 2)
          % one foot goes from in contact to not
          swing_pos = obj.footstep_contact_pos(obj.footstep_support_interval(2,:) == i);
          stance_pos = obj.footstep_contact_pos(obj.footstep_support_interval(1,:)<=i-1 & obj.footstep_support_interval(2,:)>i);
          obj.wrench_disturbance_pt(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(swing_pos{1},2) mean(swing_pos{2},2) mean(stance_pos{1},2) mean(stance_pos{2},2)],2));
          obj.com_guess(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2) mean(swing_pos{1},2) mean(swing_pos{2},2)],2)+[0;0;com_height_guess]);
        elseif(length(bodies0) == 4 && length(bodies1) == 3)
          % one foot goes from in contact to not
          swing_pos = obj.footstep_contact_pos{obj.footstep_support_interval(2,:) == i};
          stance_pos = obj.footstep_contact_pos(obj.footstep_support_interval(1,:)<=i-1 & obj.footstep_support_interval(2,:)>i);
          obj.wrench_disturbance_pt(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(swing_pos,2) mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2)],2));
          obj.com_guess(:,t_idx) = bsxfun(@times,ones(1,length(t_idx)),mean([mean(stance_pos{1},2) mean(stance_pos{2},2) mean(stance_pos{3},2) mean(swing_pos,2)],2)+[0;0;com_height_guess]);
        end
        % compute the contact wrench set
        footstep_id = find(obj.footstep_support_interval(1,:)<=i-1 & obj.footstep_support_interval(2,:)>i-1);
        footstep_cws_ray = [];
        for j = 1:length(footstep_id)
          for k = 1:size(obj.footstep_contact_pos{footstep_id(j)},2)
            footstep_cws_ray = [footstep_cws_ray [obj.footstep_fc_edges{footstep_id(j)};crossSkewSymMat(obj.footstep_contact_pos{footstep_id(j)}(:,k))*obj.footstep_fc_edges{footstep_id(j)}]];
          end
        end
        % transforms from rays to the facets representation.
        [Ain,bin] = raysToLincon(footstep_cws_ray);
        for j = 1:length(t_idx)
          obj.Ain_cws{t_idx(j)} = Ain;
          obj.bin_cws{t_idx(j)} = bin;
        end
      end
    end
  end
end