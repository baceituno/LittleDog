classdef Quad_FootstepPlan
% A container for footstep plans. A footstep plan contains the location of each
% footstep, the safe terrain regions that those steps occupy, and the assignments
% of each step to a safe region.
  properties
    footsteps % a list of Footstep objects
    params % footstep plan params, as in drc.footstep_plan_params_t
    safe_regions % a list of safe regions, as described in planFootsteps.m
    region_order % a list of the same length as footsteps. If region_order(i) == j, then footsteps(i).pos must be within safe_regions(j)
    quadruped
  end

  methods
    function obj = Quad_FootstepPlan(footsteps, quadruped, params, safe_regions, region_order)
      obj.footsteps = footsteps;
      obj.quadruped = quadruped;
      obj.params = struct(params);
      obj.safe_regions = safe_regions;
      obj.region_order = region_order;
    end

    function plan = slice(obj, idx)
      plan = obj;
      plan.footsteps = obj.footsteps(idx);
      plan.region_order = obj.region_order(idx);
    end

    function plan = concat(obj, plan_extra)
      plan = obj;
      n_steps = length(plan_extra.footsteps);
      for j=5:n_steps
        plan.footsteps(end+1) = plan_extra.footsteps(j);
        plan.region_order(end+1) = plan_extra.region_order(j);
      end
    end

    function plan = extend(obj, final_length, n)
      % Extend a footstep plan by replicating its final n footsteps. Useful for
      % generating seeds for later optimization.
      % @param final_length desired number of footsteps in the extended plan
      % @option n how many final steps to consider (the last n steps will be
      %          repeatedly appended to the footstep plan until the final
      %          length is achieved). Optional. Default: 6
      % @retval plan the extended plan
      if nargin < 3
        n = 6;
      end
      if n > length(obj.footsteps)
        error('DRC:FootstepPlan:NotEnoughStepsToExtend', 'Not enough steps in the plan to extend in the requested manner');
      end
      if final_length <= length(obj.footsteps)
        plan = plan.slice(1:final_length);
      else
        plan = obj;
        j = 1;
        source_ndx = (length(obj.footsteps) - n) + (1:n);
        for k = (length(obj.footsteps) + 1):final_length
          if mod(k,4) == 0
            j = 4;  
          elseif mod(k,3) == 0
            j = 3;  
          elseif mod(k,2) == 0
            j = 2;
          else
            j = 1;
          end
          plan.footsteps(k) = plan.footsteps(source_ndx(j));
          plan.region_order(k) = plan.region_order(source_ndx(j));
          plan.footsteps(k).id = plan.footsteps(k-1).id + 1;
          j = mod(j, length(source_ndx)) + 1;
        end
      end
    end

    function ts = compute_step_timing(obj, quadruped)
      % Compute the approximate step timing based on the distance each swing foot must travel.
      % @retval ts a vector of times (in seconds) corresponding to the completion
      %            (return to double support) of each step in the plan. The first
      %            two entries of ts will always be zero, since the first two steps
      %            in the plan correspond to the current locations of the feet.
      if nargin < 2
        quadruped = obj.quadruped;
      end
      ts = zeros(1, length(obj.footsteps));
      for j = 5:length(obj.footsteps)
        if mod(j,2)
          n = 3;
        else
          n = 1;
        end
        [swing_ts, ~, ~, ~] = planSwing(quadruped, obj.footsteps(j-n), obj.footsteps(j));
        ts(j) = ts(j-n) + swing_ts(end);
      end
    end

    function varargout = sanity_check(obj)
      % Self-test for footstep plans.
      ok = true;
      frame_ids = [obj.footsteps.frame_id];
      if any(frame_ids(1:end-1) == frame_ids(2:end))
        ok = false;
        if nargout < 1
          error('Body indices should not repeat.');
        end
      end
      varargout = {ok};
    end

    function draw_lcmgl(obj, lcmgl)
      lcmgl.glColor3f(1,1,0);
      k = 1;
      for j = 1:length(obj.footsteps)
        if mod(k, 4) == 0
          lcmgl.glColor3f(0,1,0);
          k = 0;
        elseif mod(k, 3) == 0
          lcmgl.glColor3f(0,1,1);
        elseif mod(k, 2) == 0
          lcmgl.glColor3f(1,0,0);
        else
          lcmgl.glColor3f(1,0,1);  
        end
        lcmgl.glPushMatrix();
        lcmgl.glTranslated(obj.footsteps(j).pos(1),...
                           obj.footsteps(j).pos(2),...
                           obj.footsteps(j).pos(3));
        axis = rpy2axis(obj.footsteps(j).pos(4:6));
        axis = axis([4,1:3]); % LCMGL wants [angle; axis], not [axis; angle]
        axis(1) = 180 / pi * axis(1);
        lcmgl.glRotated(axis(1), axis(2), axis(3), axis(4));
        lcmgl.sphere([0;0;0], 0.015, 20, 20);
        lcmgl.glPushMatrix();
        len = 0.05;
        lcmgl.glTranslated(len / 2, 0, 0);
        lcmgl.drawArrow3d(len, 0.02, 0.02, 0.005);
        lcmgl.glPopMatrix();
        lcmgl.glPopMatrix();
        k = k + 1;
      end
    end

    function draw_2d(obj)
      hold on
      X1 = [obj.footsteps(1:6:end).pos];
      X2 = [obj.footsteps(2:6:end).pos];
      X3 = [obj.footsteps(3:6:end).pos];
      X4 = [obj.footsteps(4:6:end).pos];
      plot(X1(1,:), X1(2,:), 'mo',...
          X2(1,:), X2(2,:), 'ro',...
          X3(1,:), X3(2,:), 'co',...
          X4(1,:), X4(2,:), 'go')
      axis equal
    end

    function steps_rel = relative_step_offsets(obj)
      % Compute the relative displacements of the footsteps (for checking collocation results from footstepNLP)
      steps = obj.step_matrix();
      nsteps = length(obj.footsteps);
      steps_rel = zeros(6, nsteps-5);
      for j = 4:nsteps
        if mod(j,2)
          n = 3;
        else
          n = 1;
        end
        R = rotmat(-steps(6,j-n));
        steps_rel(:,j-n) = [R * (steps(1:2,j) - steps(1:2,j-n));
                    steps(3:6,j) - steps(3:6,j-n)];
      end
    end

    function steps = step_matrix(obj)
      % Return the footstep positions as a 6xnsteps matrix in x y z roll
      % pitch yaw
      steps = [obj.footsteps.pos];
    end

    function obj = trim_duplicates(obj, trim_threshold)
      % Remove steps from the beginning which are identical to the initial
      % poses and remove steps at the end which are identical to the final poses,
      % returning a modified copy of the plan.
      if nargin < 2
        trim_threshold = [0.02; 0.02; 0.02; inf; inf; pi/32];
      end
      poses = [obj.footsteps.pos];
      nsteps = length(obj.footsteps);

      trim_init = false(1, nsteps);
      trim_init(1:6) = 1;
      trim_final = false(1, nsteps);
      trim_final(end-5:end) = 1;
      % initial trim
      k = 0;
      for j = 5:nsteps
        if sum(~trim_init) + 4 < obj.params.min_num_steps
          break
        end
        if mod(k, 4) == 0
          trim_init(j) = all(abs(poses(:,j) - poses(:,4)) <= trim_threshold);
        elseif mod(k, 3) == 0
          trim_init(j) = all(abs(poses(:,j) - poses(:,3)) <= trim_threshold);
        elseif mod(k, 2) == 0
          trim_init(j) = all(abs(poses(:,j) - poses(:,2)) <= trim_threshold);
        else
          trim_init(j) = all(abs(poses(:,j) - poses(:,1)) <= trim_threshold);
        end
        k = k + 1;
      end
      trim_init(find(trim_init, 4, 'last')) = false;

      % final trim
      k = 4;
      for j = 3:nsteps-4
        if sum(~trim_final) + 4 < obj.params.min_num_steps
          break
        end
        if mod(k, 4) == 0
          trim_final(j) = all(abs(poses(:,j) - poses(:,end-2)) <= trim_threshold);
        elseif mod(k, 3) == 0
          trim_final(j) = all(abs(poses(:,j) - poses(:,end-3)) <= trim_threshold);
        elseif mod(k, 2) == 0
          trim_final(j) = all(abs(poses(:,j) - poses(:,end-4)) <= trim_threshold);
        else
          trim_final(j) = all(abs(poses(:,j) - poses(:,end-5)) <= trim_threshold);
          k = 5;
        end
      end

      trim_final(find(trim_final, 4, 'first')) = false;

      obj = obj.slice(~(trim_init | trim_final));
      k = k - 1;

    end
  end
  methods(Static=true)

    function plan = blank_plan(quadruped, nsteps, ordered_frame_id, params, safe_regions)
      % Construct a FootstepPlan with all footstep poses set to NaN, but with the individual
      % step frame_id fields filled out appropriately. This is useful because all of the
      % existing footstep optimizations require that the maximum number of footsteps and the
      % sequence of the feet be assigned beforehand.
      footsteps = Quad_Footstep.empty();
      n = 1;
      for j = 1:nsteps
        pos = nan(6,1);
        id = j;
        if mod(n, 4) == 0
          n = 0;
          k = 4;
        elseif mod(n, 3) == 0
          k = 3;
        elseif mod(n, 2) == 0
          k = 2;
        else
          k = 1;
        end
        frame_id = ordered_frame_id(k);
        is_in_contact = true;
        pos_fixed = zeros(6,1);
        terrain_pts = [];
        infeasibility = nan;
        walking_params = [];
        footsteps(j) = Quad_Footstep(pos, id, frame_id, is_in_contact, pos_fixed, terrain_pts, infeasibility, walking_params);
        n = n + 1;
      end
      region_order = nan(1, nsteps);
      plan = Quad_FootstepPlan(footsteps, quadruped, params, safe_regions, region_order);
    end
  end
end
