classdef Quad_MixedIntegerFootstepPlanningProblem < Quad_MixedIntegerConvexProgram
% Structure for Quadruped footstep planning through Mixed Integer Convex Optimization
% Currently implemented with numerical variables for constraints in reachability and 
% robot shape.
%
% Developed by Bernardo Aceituno-C (Mechatronics Group, USB C Laboratory).
%
% As presented in the 2017 paper "A Generalized Mixed-Integer Convex Model
% For Multilegged Footstep Planning on Uneven Terrain" by Bernardo Aceituno-Cabezas,
% J. Cappelletto, J. Grieco, and G. Fernarndez (USB)
%
% Based on MixedIntegerFootstepPlanningProblem by Robin Deits (Robot Locomotion Group, MIT CSAIL).

  properties
    Quadruped;
    n_legs = 4;
    nsteps;
    seed_plan;
    plan_gait = false;
    L_leg = 0.0991;
    offset = [0.8339, -0.8339, pi - 0.8339, -pi + 0.8339];
    d_lim = 0.01;
    l_bnd = 0.1;
    gamma_z = 0.005;
    dz = 0.05;
    weights;
    max_distance = 500;
    pose_indices = [1,2,3,6];
  end

  methods
      function obj = Quad_MixedIntegerFootstepPlanningProblem(Quadruped, seed_plan)
        % Constructs the optimization problem and declares the variables for each footstep
        % @param Quadruped a Quadruped.
        % @param seed_plan a blank footstep plan, provinding the structure of the
        %                  desired plan. Probably generated with
        %                  FootstepPlan.blank_plan()

        typecheck(Quadruped, 'Quadruped');
        typecheck(seed_plan, 'Quad_FootstepPlan');

        obj = obj@Quad_MixedIntegerConvexProgram();

        obj.Quadruped = Quadruped;
        obj.seed_plan = seed_plan;
        obj.nsteps = length(obj.seed_plan.footsteps);
        obj.n_legs = 4;

        assert(mod(obj.nsteps,obj.n_legs) == 0)
        
        obj.weights = obj.Quadruped.getFootstepOptimizationWeights();

        seed_steps = [seed_plan.footsteps.pos];
        min_yaw = -pi;
        max_yaw = pi;

        lb = [repmat(seed_steps(1:3,1) - obj.max_distance, 1, obj.nsteps);
              min_yaw + zeros(1, obj.nsteps)];
        ub = [repmat(seed_steps(1:3,1) + obj.max_distance, 1, obj.nsteps);
              max_yaw + zeros(1, obj.nsteps)];
        
        lb(:,1) = seed_steps(obj.pose_indices, 1);
        lb(:,2) = seed_steps(obj.pose_indices, 2);
        lb(:,3) = seed_steps(obj.pose_indices, 3);
        lb(:,4) = seed_steps(obj.pose_indices, 4);

        ub(:,1) = seed_steps(obj.pose_indices, 1);
        ub(:,2) = seed_steps(obj.pose_indices, 2);
        ub(:,3) = seed_steps(obj.pose_indices, 3);
        ub(:,4) = seed_steps(obj.pose_indices, 4);
        
        obj = obj.addVariable('footsteps', 'C', [obj.n_legs, obj.nsteps], lb, ub);
      end

      function obj = addQuadraticGoalObjective(obj, goal_pose, step_indices)
        % For each index j in step_indices, add a caudratic cost of the form:
        % weights(j) * (footstep(j) - xgoal(6 - N +j))' * Qfoal * (footstep(j) - xgoal(6 - N +j)
        
        k = 1;
        w_goal = diag(obj.weights.goal(obj.pose_indices));
        for j = step_indices
          if mod(k,4) == 0
            xg = goal_pose.Foot4(obj.pose_indices);
          elseif mod(k,3) == 0
            xg = goal_pose.Foot3(obj.pose_indices);
          elseif mod(k,2) == 0
            xg = goal_pose.Foot2(obj.pose_indices);
          else
            xg = goal_pose.Foot1(obj.pose_indices);
          end
          Qi = sparse([], [], [], obj.nv, obj.nv, 4);
          ci = zeros(obj.nv, 1);
          Qi(obj.vars.footsteps.i(:,j), obj.vars.footsteps.i(:,j)) = w_goal;
          ci(obj.vars.footsteps.i(:,j)) = -2 * w_goal * xg;
          objcon_i = xg' * w_goal * xg;
          obj = obj.addCost(Qi, ci, objcon_i);
          k = k + 1;
        end
      end

      function obj = addTrimToFinalPoses(obj)
        % Add a binary variable for each footstep which, if true, forces that footstep to the
        % final pose in the footstep plan. This allows us to trim it out of the footstep plan later.
        % A linear objective placed on those trim variables lets us tune the number of footsteps
        % in the plan.
        
        obj = obj.addVariable('trim', 'B', [1, obj.nsteps], 0, 1);
        w_trim = 1;
        min_num_steps = max([obj.seed_plan.params.min_num_steps + 6, 7]);

        obj.vars.trim.lb(end-3:end) = 1;
        obj.vars.trim.ub(end-3:end) = 1;
        obj.vars.trim.lb(1:4) = 0;
        obj.vars.trim.ub(1:4) = 0;

        Ai = sparse(obj.nsteps-4 + 4 + max(obj.nsteps-8, 0) * 8, obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset = 0;
        expected_offset = size(Ai, 1);
        for j = 2:obj.nsteps
          % trim(j) >= trim(j-1)
          Ai(offset+1, obj.vars.trim.i(j)) = -1;
          Ai(offset+1, obj.vars.trim.i(j-1)) = 1;
          offset = offset + 1;
        end
        % sum(trim) <= obj.nsteps - (min_num_steps - 4)
        Ai(offset+1, obj.vars.trim.i) = 1;
        bi(offset+1) = obj.nsteps - (min_num_steps - 4);
        offset = offset + 1;
        M = obj.max_distance;
        
        for j = 5:obj.nsteps
          % x(:,j) - x(:,k) <= M(1-trim(j))
          Ai(offset+(1:4), obj.vars.footsteps.i(:,j)) = speye(4);
          Ai(offset+(1:4), obj.vars.footsteps.i(:,j-4)) = -speye(4);
          Ai(offset+(1:4), obj.vars.trim.i(j)) = M;
          bi(offset+(1:4)) = M;
          offset = offset + 4;

          % x(:,j) - x(:,k) >= -M(1-trim(j))
          Ai(offset+(1:4), obj.vars.footsteps.i(:,j)) = -speye(4);
          Ai(offset+(1:4), obj.vars.footsteps.i(:,j-4)) = speye(4);
          Ai(offset+(1:4), obj.vars.trim.i(j)) = M;
          bi(offset+(1:4)) = M;
          offset = offset + 4; 
        end

        obj = obj.addLinearConstraints(Ai, bi, [], []);
        assert(offset == expected_offset);

        c = zeros(obj.nv, 1);
        c(obj.vars.trim.i) = -w_trim;
        obj = obj.addCost([], c, w_trim * obj.nsteps);
      end

      function obj = addQuadraticRelativeObjective(obj)
        % Add a quadratic cost on the transfer distance between configurations
        % in the form:
        % (pi - pi-n)'Qr(pi - pi-n)
        % Where Qr is chosen to matain the footsteps displacement on its nominal position

        % defines the center of contacts as the geometric center of each n_leg configuration

        obj = obj.addVariable('p', 'C', [3, obj.nsteps], -10*ones(3,obj.nsteps), 10*ones(3 ,obj.nsteps));
        k = 1;
        for j = 1:obj.nsteps
          if mod(k,4) == 0
            k = 0;
            n1 = j-3;
            n2 = j-2;
            n3 = j-1;
            n4 = j;
            correct = obj.offset(4);
          elseif mod(k,3) == 0
            n1 = j-2;
            n2 = j-1;
            n3 = j;
            n4 = j+1;
            correct = obj.offset(3);
          elseif mod(k,2) == 0
            n1 = j-1;
            n2 = j;
            n3 = j+1;
            n4 = j+2;
            correct = obj.offset(2);
          else
            n1 = j;
            n2 = j+1;
            n3 = j+2;
            n4 = j+3;
            correct = obj.offset(1);
          end
          %setting the values
          Aeq = sparse(3,obj.nv);
          beq = zeros(3,1);

          Aeq(1,obj.vars.footsteps.i(1,n1)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n2)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n3)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n4)) = 1;
          Aeq(1,obj.vars.p.i(1,j)) = -4;

          Aeq(2,obj.vars.footsteps.i(2,n1)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n2)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n3)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n4)) = 1;
          Aeq(2,obj.vars.p.i(2,j)) = -4;

          Aeq(3,obj.vars.footsteps.i(3,n1)) = 1;
          Aeq(3,obj.vars.footsteps.i(3,n2)) = 1;
          Aeq(3,obj.vars.footsteps.i(3,n3)) = 1;
          Aeq(3,obj.vars.footsteps.i(3,n4)) = 1;
          Aeq(3,obj.vars.p.i(3,j)) = -4;
            
          k = k + 1;

          obj = obj.addLinearConstraints([],[],Aeq,beq);
        end
        
        % adds a quadratic cost in the relative distance between configurations
        % defined as the transfer distance or T_d(j) = p(j) - p(j-n_legs)
        % and a quadratic cost in the distance between each step and the center of each safe region

        Qr = 1;
        
        for j = 5:obj.nsteps
          Qi = sparse(obj.nv, obj.nv);

          Qi(obj.vars.p.i(1,j), obj.vars.p.i(1,j)) = Qr;
          Qi(obj.vars.p.i(1,j), obj.vars.p.i(1,j-4)) = -Qr;
          Qi(obj.vars.p.i(1,j-4), obj.vars.p.i(1,j)) = -Qr;
          Qi(obj.vars.p.i(1,j-4), obj.vars.p.i(1,j-4)) = Qr;

          obj = obj.addCost(Qi, [], []);
          Qi = sparse(obj.nv, obj.nv);

          Qi(obj.vars.p.i(2,j), obj.vars.p.i(2,j)) = Qr;
          Qi(obj.vars.p.i(2,j), obj.vars.p.i(2,j-4)) = -Qr;
          Qi(obj.vars.p.i(2,j-4), obj.vars.p.i(2,j)) = -Qr;
          Qi(obj.vars.p.i(2,j-4), obj.vars.p.i(2,j-4)) = Qr;

          obj = obj.addCost(Qi, [], []);
          Qi = sparse(obj.nv, obj.nv);

          Qi(obj.vars.footsteps.i(1,j), obj.vars.footsteps.i(1,j)) = Qr;
          Qi(obj.vars.footsteps.i(1,j), obj.vars.r_nom.i(1,j)) = -Qr;
          Qi(obj.vars.r_nom.i(1,j), obj.vars.footsteps.i(1,j)) = -Qr;
          Qi(obj.vars.r_nom.i(1,j), obj.vars.r_nom.i(1,j)) = Qr;

          obj = obj.addCost(Qi, [], []);
          Qi = sparse(obj.nv, obj.nv);

          Qi(obj.vars.footsteps.i(2,j), obj.vars.footsteps.i(2,j)) = Qr;
          Qi(obj.vars.footsteps.i(2,j), obj.vars.r_nom.i(2,j-4)) = -Qr;
          Qi(obj.vars.r_nom.i(2,j-4), obj.vars.footsteps.i(2,j)) = -Qr;
          Qi(obj.vars.r_nom.i(2,j-4), obj.vars.r_nom.i(2,j-4)) = Qr;

          obj = obj.addCost(Qi, [], []);
          Qi = sparse(obj.nv, obj.nv);

          Qi(obj.vars.footsteps.i(4,j), obj.vars.footsteps.i(4,j)) = Qr;
          Qi(obj.vars.footsteps.i(4,j), obj.vars.footsteps.i(4,j-4)) = -Qr;
          Qi(obj.vars.footsteps.i(4,j-4), obj.vars.footsteps.i(4,j)) = -Qr;
          Qi(obj.vars.footsteps.i(4,j-4), obj.vars.footsteps.i(4,j-4)) = Qr;

          obj = obj.addCost(Qi, [], []);
        end
      end

      function obj = addGaitConstraints(obj)
        % Add mixed-integer constraints that require that assign each
        % leg movement to a timeslot between 1 and nsteps
        % such that for footstep i T(i,j) implies that leg i moves in timestep j
        % where T is a binary matrix with sum(T(i)) == 1

        obj.plan_gait = true;

        obj = obj.addVariable('seq', 'B', [obj.nsteps-4, obj.nsteps-4], 0, 1);
        obj = obj.addVariable('timming', 'C', [1, obj.nsteps-4], 1, obj.nsteps);

        obj.vars.seq.lb(1,1) = 1;
        obj.vars.seq.ub(1,1) = 1;

        for i = 1:obj.nsteps-4
          % assigns each leg movement to a time variable
          Aeq = sparse(1,obj.nv);
          beq = 0;
          Aeq(1,obj.vars.seq.i(i,:)) = 1:obj.nsteps-4;
          Aeq(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          % constrains that that each step is taken in only
          % one timminglot
          Aeq = sparse(1,obj.nv);
          beq = 1;
          Aeq(1,obj.vars.seq.i(i,:)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end

        % constrains that each leg moves after the previous configuration has
        % completely transfered to the current configuration
        k = 1;
        for i = 5:obj.nsteps-4
          if k > 4
            k = 1;
          end
          Ai = sparse(1,obj.nv);
          bi = -1;
          Ai(1,obj.vars.timming.i(1,i-3-k)) = 1;
          Ai(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(1,obj.nv);
          bi = -1;
          Ai(1,obj.vars.timming.i(1,i-2-k)) = 1;
          Ai(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(1,obj.nv);
          bi = -1;
          Ai(1,obj.vars.timming.i(1,i-1-k)) = 1;
          Ai(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(1,obj.nv);
          bi = -1;
          Ai(1,obj.vars.timming.i(1,i-k)) = 1;
          Ai(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(1,obj.nv);
          bi = -2;
          Ai(1,obj.vars.timming.i(1,i-4)) = 1;
          Ai(1,obj.vars.timming.i(1,i)) = -1;
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          k = k + 1;
        end

        % constrains the possible configurations in the gait
        for i = 1:4:obj.nsteps-7
          for j = 1:obj.nsteps-4
            Ai = sparse(1,obj.nv);
            bi = 1;
            Ai(1,obj.vars.seq.i(i,j)) = 1;
            Ai(1,obj.vars.seq.i(i+1,j)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Ai = sparse(1,obj.nv);
            Ai(1,obj.vars.seq.i(i+2,j)) = 1;
            Ai(1,obj.vars.seq.i(i+3,j)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Ai = sparse(1,obj.nv);
            Ai(1,obj.vars.seq.i(i,j)) = 1;
            Ai(1,obj.vars.seq.i(i+2,j)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);

            Ai = sparse(1,obj.nv);
            Ai(1,obj.vars.seq.i(i+1,j)) = 1;
            Ai(1,obj.vars.seq.i(i+3,j)) = 1;
            obj = obj.addLinearConstraints(Ai, bi, [], []);
          end
        end

        % defines variables for the height deviation on each configuration
        % and for the binary approximation of the variable
        gamma_z = obj.gamma_z;
        M = 100; % Large height difference for big-M formulation

        obj = obj.addVariable('dh', 'B', [2, obj.nsteps-4], 0, 1);
        obj = obj.addVariable('desv', 'C', [1, obj.nsteps-4], -inf, inf);
        obj = obj.addVariable('delta_z', 'C', [1, obj.nsteps-4], -inf, inf);

        for j = 5:obj.nsteps
          % computes the height deviation for each configuration
          Aeq = sparse(1,obj.nv);
          beq = 0;
          Aeq(1,obj.vars.desv.i(1,j-4)) = -1;
          Aeq(1,obj.vars.footsteps.i(3,j)) = 10;
          Aeq(1,obj.vars.footsteps.i(3,j-4)) = -10;
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          % makes a dissagregate approximation of the deviation magnitude
          Ai = sparse(2,obj.nv);
          bi = [M - gamma_z; gamma_z];

          Ai(1, obj.vars.desv.i(1,j-4)) = -1;
          Ai(2, obj.vars.desv.i(1,j-4)) = 1;
          
          Ai(1, obj.vars.dh.i(1,j-4)) = M;
          Ai(2, obj.vars.dh.i(1,j-4)) = -M;
          
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Ai = sparse(2,obj.nv);
          bi = [M - gamma_z; gamma_z];

          Ai(1, obj.vars.desv.i(1,j-4)) = 1;
          Ai(2, obj.vars.desv.i(1,j-4)) = -1;
          
          Ai(1, obj.vars.dh.i(2,j-4)) = M;
          Ai(2, obj.vars.dh.i(2,j-4)) = -M;
          
          obj = obj.addLinearConstraints(Ai, bi, [], []);

          Aeq = sparse(1,obj.nv);
          beq = 0;
          Aeq(1, obj.vars.dh.i(1,j-4)) = 1;
          Aeq(1, obj.vars.dh.i(2,j-4)) = 1;
          Aeq(1, obj.vars.delta_z.i(1,j-4)) = -1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);
        end

        obj = obj.addVariable('s', 'C', [obj.nsteps-4, obj.nsteps-4], -inf, inf);
        obj = obj.addVariable('cost_desv', 'C', [1, obj.nsteps-4], -inf, inf);

        % adds a cost on the duration of the plan
        q_t = 2;
        ci = zeros(obj.nv, 1);
        for j = 1:obj.nsteps-4
          ci(obj.vars.timming.i(1,j),1) = q_t;
        end
        obj = obj.addCost([], ci, []);

        q_d = 100;
        M = 1000000;

        % adds a cuadratic cost on the sum of legs in transfer on
        % the unevenness of the terrain
        for j = 1:obj.nsteps-4
          Qi = sparse(obj.nv, obj.nv);
          ci = sparse(obj.nv, 1);

          for i = 1:obj.nsteps-4
            % slack_ij = seq(i,j)*dh(i)
            Ain = sparse(4, obj.nv);
            bin = zeros(4,1);
            Ain(1,obj.vars.seq.i(i,j)) = -M;
            Ain(1,obj.vars.s.i(i,j)) = 1;
            bin(1,1) = 0;

            Ain(2,obj.vars.seq.i(i,j)) = -M;
            Ain(2,obj.vars.s.i(i,j)) = -1;
            bin(2,1) = 0;

            Ain(3,obj.vars.seq.i(i,j)) = M;
            Ain(3,obj.vars.s.i(i,j)) = 1;
            Ain(3,obj.vars.delta_z.i(1,i)) = -1;
            bin(3,1) = M;

            Ain(4,obj.vars.seq.i(i,j)) = M;
            Ain(4,obj.vars.s.i(i,j)) = -1;
            Ain(4,obj.vars.delta_z.i(1,i)) = 1;
            bin(4,1) = M;

            obj = obj.addLinearConstraints(Ain, bin, [], []);
          end

          Aeq = sparse(1, obj.nv);
          beq = zeros(1,1);

          Aeq(1,obj.vars.cost_desv.i(1,j)) = -1;
          Aeq(1,obj.vars.s.i(:,j)) = 1;
          obj = obj.addLinearConstraints([], [], Aeq, beq);

          Qi(obj.vars.cost_desv.i(1,j), obj.vars.cost_desv.i(1,j)) = q_d;
          obj = obj.addCost(Qi, [], []);
        end
      end

      function obj = addTerrainRegions(obj, safe_regions)
        % Add mixed-integer constraints that require that 
        % each footstep lie within one of those safe regions described by A and b.
        % such that for footstep i H(i,:) implies that A*fi < b
        % where H is a binary matrix with sum(H(i)) == 1
        
        if isempty(safe_regions)
          safe_regions = obj.seed_plan.safe_regions;
        end
        nr = length(safe_regions);
        obj = obj.addVariable('region', 'B', [nr, obj.nsteps], 0, 1);
        obj = obj.addVariable('cent_region', 'C', [3, obj.nsteps], 0, 1);

        Ai = sparse((obj.nsteps-4) * sum(cellfun(@(x) size(x, 1) + 2, {safe_regions.A})), obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset_ineq = 0;
        Aeq = sparse(obj.nsteps-4, obj.nv);
        beq = ones(obj.nsteps-4, 1);
        offset_eq = 0;

        for r = 1:nr
          A = safe_regions(r).A;
          b = safe_regions(r).b;
          Ar = [A(:,1:2), sparse(size(A, 1), 1), 0*A(:,3);
                safe_regions(r).normal', 0;
                -safe_regions(r).normal', 0];
          br = [b;
                safe_regions(r).normal' * safe_regions(r).point;
                -safe_regions(r).normal' * safe_regions(r).point];
          s = size(Ar, 1);
          M = obj.max_distance;
          for j = 5:obj.nsteps
            Ai(offset_ineq + (1:s), obj.vars.footsteps.i(:,j)) = Ar;
            Ai(offset_ineq + (1:s), obj.vars.region.i(r,j)) = M;
            bi(offset_ineq + (1:s)) = br + M;
            offset_ineq = offset_ineq + s;
          end
        end
        assert(offset_ineq == size(Ai, 1));
        for j = 5:obj.nsteps
          Aeq(offset_eq + 1, obj.vars.region.i(:,j)) = 1;
          offset_eq = offset_eq + 1;
        end
        assert(offset_eq == size(Aeq, 1));
        obj = obj.addLinearConstraints(Ai, bi, Aeq, beq);
      end

      function obj = fixRotation(obj)
        % Fix the rotations of every step to the initial value (no rotations allowed)
        k = 1;
        yaw = obj.seed_plan.footsteps(1).pos(6);
        for j = 5:obj.nsteps
          obj.seed_plan.footsteps(j).pos(6) = yaw;
          obj.vars.footsteps.lb(4,j) = yaw - 0.1;
          obj.vars.footsteps.ub(4,j) = yaw + 0.1;
          k = k + 1;
        end

        seed_steps = [obj.seed_plan.footsteps.pos];
        yaw = seed_steps(6,:);
        obj = obj.addVariable('cos_yaw', 'C', [1, obj.nsteps], cos(yaw), cos(yaw));
        obj = obj.addVariable('sin_yaw', 'C', [1, obj.nsteps], sin(yaw), sin(yaw));
      end

      function obj = addReachabilityConstraints(obj)
        % Adds a set of linear constraints on the reachability of each leg and the 
        % geometry of the robot.
        %
        % as presented in the 2017 paper "A Generalized Mixed-Integer Convex Model
        % For Multilegged Footstep Planning on Uneven Terrain" by Bernardo Aceituno-Cabezas,
        % J. Cappelletto, J. Grieco, and G. Fernarndez (USB)

        % Adds a linear a constraint to matain each footstep in a squared region centered on its reference 
        % feet position with respect to the Center of Contacts (CoC) with a side d_lim.

        % defines the variable r_nom for the reference position of each step
        % measured from the center of contacts of each configuration
        obj = obj.addVariable('r_nom', 'C', [2, obj.nsteps], -10*ones(2,obj.nsteps), 10*ones(2,obj.nsteps));
        
        % defines the conter of contacts as p(i) = [sum(f(1,i));sum(f(2,i))]/4
        % and the reference footstep as:
        %  r_nom(i,1) = p(i,1) + L*cos(yaw(i))
        %  r_nom(i,2) = p(i,2) + L*sin(yaw(i))

        k = 1;
        for j = 1:obj.nsteps
          if mod(k,4) == 0
            k = 0;
            n1 = j-3;
            n2 = j-2;
            n3 = j-1;
            n4 = j;
            correct = obj.offset(4);
          elseif mod(k,3) == 0
            n1 = j-2;
            n2 = j-1;
            n3 = j;
            n4 = j+1;
            correct = obj.offset(3);
          elseif mod(k,2) == 0
            n1 = j-1;
            n2 = j;
            n3 = j+1;
            n4 = j+2;
            correct = obj.offset(2);
          else
            n1 = j;
            n2 = j+1;
            n3 = j+2;
            n4 = j+3;
            correct = obj.offset(1);
          end

          %setting the values
          Aeq = sparse(2,obj.nv);
          beq = zeros(2,1);

          Aeq(1,obj.vars.footsteps.i(1,n1)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n2)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n3)) = 1;
          Aeq(1,obj.vars.footsteps.i(1,n4)) = 1;
          Aeq(1,obj.vars.r_nom.i(1,j)) = -4;
          Aeq(1,obj.vars.cos_yaw.i(n1)) = 4*obj.L_leg*cos(correct);
          Aeq(1,obj.vars.sin_yaw.i(n1)) = -4*obj.L_leg*sin(correct);

          Aeq(2,obj.vars.footsteps.i(2,n1)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n2)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n3)) = 1;
          Aeq(2,obj.vars.footsteps.i(2,n4)) = 1;
          Aeq(2,obj.vars.r_nom.i(2,j)) = -4;
          Aeq(2,obj.vars.cos_yaw.i(n1)) = 4*obj.L_leg*sin(correct);
          Aeq(2,obj.vars.sin_yaw.i(n1)) = 4*obj.L_leg*cos(correct);

          k = k + 1;

          obj = obj.addLinearConstraints([],[],Aeq,beq);
        end

        k = 1;

        %we now have that f(i) must be in a square of side d_lim from the step
        
        %polytopic relaxation of the geometric contraint
        for j = 5:obj.nsteps
          Ai = sparse(4,obj.nv);
          bi = zeros(4,1);

          Ai(1,obj.vars.footsteps.i(1,j)) = 1;
          Ai(1,obj.vars.r_nom.i(1,j)) = -1;
          bi(1,1) = obj.d_lim;

          Ai(2,obj.vars.footsteps.i(1,j)) = -1;
          Ai(2,obj.vars.r_nom.i(1,j)) = 1;
          bi(2,1) = obj.d_lim;

          Ai(3,obj.vars.footsteps.i(2,j)) = 1;
          Ai(3,obj.vars.r_nom.i(2,j)) = -1;
          bi(3,1) = obj.d_lim;

          Ai(4,obj.vars.footsteps.i(2,j)) = -1;
          Ai(4,obj.vars.r_nom.i(2,j)) = 1;
          bi(4,1) = obj.d_lim;

          obj = obj.addLinearConstraints(Ai,bi,[],[]);
        end

        % Constraints the reachability of each leg for the robot configuration and its successor
        % Adds linear constraints such that every footstep lies on a square of side l_bnd 
        % centered on the previous nominal leg position.
        
        for j = 5:obj.nsteps
          %defines the linear constraint on reachability on x
          Ai = sparse(2,obj.nv);
          bi = zeros(2,1);

          Ai(1,obj.vars.footsteps.i(1,j)) = 1;
          Ai(1,obj.vars.r_nom.i(1,j-4)) = -1;
          bi(1,1) = obj.l_bnd;

          Ai(2,obj.vars.footsteps.i(1,j)) = -1;
          Ai(2,obj.vars.r_nom.i(1,j-4)) = 1;
          bi(2,1) = obj.l_bnd;

          %adds the constraints
          obj = obj.addLinearConstraints(Ai,bi,[],[]);
          
          %defines the linear constraint on reachability on y
          Ai = sparse(2,obj.nv);
          bi = zeros(2,1);

          Ai(1,obj.vars.footsteps.i(2,j)) = 1;
          Ai(1,obj.vars.r_nom.i(2,j-4)) = -1;
          bi(1,1) = obj.l_bnd;
          
          Ai(2,obj.vars.footsteps.i(2,j)) = -1;
          Ai(2,obj.vars.r_nom.i(2,j-4)) = 1;
          bi(2,1) = obj.l_bnd;

          %adds the constraints
          obj = obj.addLinearConstraints(Ai,bi,[],[]);
          
          %defines the linear constraint on reachability of the yaw
          Ai = sparse(2,obj.nv);
          bi = zeros(2,1);

          Ai(1,obj.vars.footsteps.i(4,j)) = 1;
          Ai(1,obj.vars.footsteps.i(4,j-4)) = -1;
          bi(1,1) = pi/8;

          Ai(2,obj.vars.footsteps.i(4,j)) = -1;
          Ai(2,obj.vars.footsteps.i(4,j-4)) = 1;
          bi(2,1) = pi/8;

          %adds the constraints
          obj = obj.addLinearConstraints(Ai,bi,[],[]);
        end
        
        % Limits the difference in z between footsteps
        Ai = sparse((obj.nsteps-4)*2, obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset = 0;
        expected_offset = size(Ai, 1);
        
        for j = 5:obj.nsteps
          % obj.symbolic_constraints = [obj.symbolic_constraints,...
          %   -obj.seed_plan.params.nom_downward_step <= x(3,j) - x(3,j-1) <= obj.seed_plan.params.nom_upward_step];
          Ai(offset+1, obj.vars.footsteps.i(3,j)) = -1;
          Ai(offset+1, obj.vars.footsteps.i(3,j-4)) = 1;
          bi(offset+1) = obj.dz;
          Ai(offset+2, obj.vars.footsteps.i(3,j)) = 1;
          Ai(offset+2, obj.vars.footsteps.i(3,j-4)) = -1;
          bi(offset+2) = obj.dz;
          offset = offset + 2;
        end

        assert(offset == expected_offset);
        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      function obj = addSinCosLinearEquality(obj)
        % Adds a linear approximation of the sin and cos functions for eacg footstep yaw 
        % as presented in the 2014 paper "Footstep Planning on Uneven Terrain with
        % Mixed-Integer Convex Optimization" by Robin Deits and Russ Tedrake (MIT)

        yaw0 = obj.seed_plan.footsteps(1).pos(6);
        min_yaw = yaw0 - pi/8 * floor((obj.nsteps-4)/4);
        max_yaw = yaw0 + pi/8 * ceil((obj.nsteps-4)/4);

        if max_yaw > pi
          max_yaw = pi;
        end

        if min_yaw < pi
          min_yaw = -pi;
        end

        cos_boundaries = reshape(bsxfun(@plus, [min_yaw:pi:max_yaw; min_yaw:pi:max_yaw], [-(pi/2-1); (pi/2-1)]), 1, []);
        sin_boundaries = reshape(bsxfun(@plus, [min_yaw:pi:max_yaw; min_yaw:pi:max_yaw], [-1; 1]), 1, []);

        obj = obj.addVariableIfNotPresent('cos_yaw', 'C', [1, obj.nsteps], -1, 1);
        obj = obj.addVariableIfNotPresent('sin_yaw', 'C', [1, obj.nsteps], -1, 1);
        obj = obj.addVariable('cos_sector', 'B', [length(cos_boundaries)-1, obj.nsteps], 0, 1);
        obj = obj.addVariable('sin_sector', 'B', [length(sin_boundaries)-1, obj.nsteps], 0, 1);
        obj = obj.addInitialSinCosConstraints();

        obj.vars.footsteps.lb(4,7:end) = min_yaw;
        obj.vars.footsteps.ub(4,7:end) = max_yaw;

        for j = 5:4:obj.nsteps-4
          %fix all legs to always point in the same direction
          n1 = j;
          n2 = j+1;
          n3 = j+2;
          n4 = j+3;

          Aeq = sparse(9, obj.nv);
          beq = zeros(9, 1);

          Aeq(1,obj.vars.footsteps.i(4,n1)) = 1;
          Aeq(1,obj.vars.footsteps.i(4,n2)) = -1;

          Aeq(2,obj.vars.footsteps.i(4,n1)) = 1;
          Aeq(2,obj.vars.footsteps.i(4,n3)) = -1;

          Aeq(3,obj.vars.footsteps.i(4,n1)) = 1;
          Aeq(3,obj.vars.footsteps.i(4,n4)) = -1;

          Aeq(4,obj.vars.cos_yaw.i(n1)) = 1;
          Aeq(4,obj.vars.cos_yaw.i(n2)) = -1;

          Aeq(5,obj.vars.cos_yaw.i(n1)) = 1;
          Aeq(5,obj.vars.cos_yaw.i(n3)) = -1;

          Aeq(6,obj.vars.cos_yaw.i(n1)) = 1;
          Aeq(6,obj.vars.cos_yaw.i(n4)) = -1;

          Aeq(7,obj.vars.sin_yaw.i(n1)) = 1;
          Aeq(7,obj.vars.sin_yaw.i(n2)) = -1;

          Aeq(8,obj.vars.sin_yaw.i(n1)) = 1;
          Aeq(8,obj.vars.sin_yaw.i(n3)) = -1;

          Aeq(9,obj.vars.sin_yaw.i(n1)) = 1;
          Aeq(9,obj.vars.sin_yaw.i(n4)) = -1;

          obj = obj.addLinearConstraints([], [], Aeq, beq); 
        end

        obj = obj.addVariable('unit_circle_slack', 'C', [1,1], norm([pi/4;pi/4]), norm([pi/4;pi/4]));
        Aeq_s = sparse(obj.nsteps, obj.nv);
        Aeq_c = zeros(obj.nsteps, obj.nv);
        beq = ones(size(Aeq_s, 1), 1);
        
        for j = 1:obj.nsteps
          Aeq_c(j, obj.vars.cos_sector.i(:,j)) = 1;
          Aeq_s(j, obj.vars.sin_sector.i(:,j)) = 1;
        end

        obj = obj.addLinearConstraints([], [], [Aeq_s; Aeq_c], [beq; beq]);
        obj = obj.addPolyConesByIndex([repmat(obj.vars.unit_circle_slack.i, 1, obj.nsteps-4); obj.vars.cos_yaw.i(5:end); obj.vars.sin_yaw.i(5:end)], 8);

        M = 2*pi;
        Ai = sparse((obj.nsteps-4) * (length(cos_boundaries)-1 + length(sin_boundaries)-1) * 4, obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset = 0;
        expected_offset = size(Ai, 1);

        for s = 1:length(cos_boundaries)-1
          th0 = cos_boundaries(s);
          th1 = cos_boundaries(s+1);

          th = (th0 + th1)/2;
          cos_slope = -sin(th);
          cos_intercept = cos(th) - (cos_slope * th);

          for j = 5:obj.nsteps
            % -yaw(j) <= -th0 + M(1-cos_sector(s,j))
            Ai(offset+1, obj.vars.footsteps.i(4,j)) = -1;
            Ai(offset+1, obj.vars.cos_sector.i(s,j)) = M;
            bi(offset+1) = -th0 + M;
            % yaw(j) <= th1 + M(1-cos_sector(s,j))
            Ai(offset+2, obj.vars.footsteps.i(4,j)) = 1;
            Ai(offset+2, obj.vars.cos_sector.i(s,j)) = M;
            bi(offset+2) = th1 + M;
            offset = offset + 2;

            % cos_yaw(j) <= cos_slope * yaw(j) + cos_intercept + M(1-cos_sector(s,j))
            Ai(offset+1, obj.vars.cos_yaw.i(j)) = 1;
            Ai(offset+1, obj.vars.footsteps.i(4,j)) = -cos_slope;
            Ai(offset+1, obj.vars.cos_sector.i(s,j)) = M;
            bi(offset+1) = cos_intercept + M;
            % cos_yaw(j) >= cos_slope * yaw(j) + cos_intercept - M(1-cos_sector(s,j))
            Ai(offset+2, obj.vars.cos_yaw.i(j)) = -1;
            Ai(offset+2, obj.vars.footsteps.i(4,j)) = cos_slope;
            Ai(offset+2, obj.vars.cos_sector.i(s,j)) = M;
            bi(offset+2) = -cos_intercept + M;
            offset = offset + 2;
          end

        end

        for s = 1:length(sin_boundaries)-1
          th0 = sin_boundaries(s);
          th1 = sin_boundaries(s+1);

          th = (th0 + th1)/2;
          sin_slope = cos(th);
          sin_intercept = sin(th) - (sin_slope * th);

          for j = 5:obj.nsteps
            % -yaw(j) <= -th0 + M(1-sin_sector(s,j))
            Ai(offset+1, obj.vars.footsteps.i(4,j)) = -1;
            Ai(offset+1, obj.vars.sin_sector.i(s,j)) = M;
            bi(offset+1) = -th0 + M;
            % yaw(j) <= th1 + M(1-sin_sector(s,j))
            Ai(offset+2, obj.vars.footsteps.i(4,j)) = 1;
            Ai(offset+2, obj.vars.sin_sector.i(s,j)) = M;
            bi(offset+2) = th1 + M;
            offset = offset + 2;

            % sin_yaw(j) <= sin_slope * yaw(j) + sin_intercept + M(1-sin_sector(s,j))
            Ai(offset+1, obj.vars.sin_yaw.i(j)) = 1;
            Ai(offset+1, obj.vars.footsteps.i(4,j)) = -sin_slope;
            Ai(offset+1, obj.vars.sin_sector.i(s,j)) = M;
            bi(offset+1) = sin_intercept + M;
            % sin_yaw(j) >= sin_slope * yaw(j) + sin_intercept - M(1-sin_sector(s,j))
            Ai(offset+2, obj.vars.sin_yaw.i(j)) = -1;
            Ai(offset+2, obj.vars.footsteps.i(4,j)) = sin_slope;
            Ai(offset+2, obj.vars.sin_sector.i(s,j)) = M;
            bi(offset+2) = -sin_intercept + M;
            offset = offset + 2;
          end
        end
        assert(offset == expected_offset);
        obj = obj.addLinearConstraints(Ai, bi, [], []);

        % Consistency between sin and cos sectors
        Ai = sparse((obj.nsteps-4) * obj.vars.sin_sector.size(1) * 2, obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset = 0;
        expected_offset = size(Ai, 1);
        nsectors = obj.vars.sin_sector.size(1);
        for k = 1:nsectors
          for j = 5:obj.nsteps
            Ai(offset+1, obj.vars.cos_sector.i(k,j)) = 1;
            Ai(offset+1, obj.vars.sin_sector.i(max(1,k-1):min(k+1,nsectors),j)) = -1;
            Ai(offset+2, obj.vars.sin_sector.i(k,j)) = 1;
            Ai(offset+2, obj.vars.cos_sector.i(max(1,k-1):min(k+1,nsectors),j)) = -1;
            offset = offset + 2;
          end
        end
        assert(offset == expected_offset);
        obj = obj.addLinearConstraints(Ai, bi, [], []);

        % Transitions between sectors
        Ai = sparse((obj.nsteps-4) * (nsectors-1) * 2, obj.nv);
        bi = zeros(size(Ai, 1), 1);
        offset = 0;
        expected_offset = size(Ai, 1);
        for j = 5:obj.nsteps
          if mod(j,2)
            for k = 1:nsectors - 1
              Ai(offset+1, obj.vars.cos_sector.i(k,j-3)) = 1;
              Ai(offset+1, obj.vars.cos_sector.i(k:k+1,j)) = -1;
              Ai(offset+2, obj.vars.sin_sector.i(k,j-3)) = 1;
              Ai(offset+2, obj.vars.sin_sector.i(k:k+1,j)) = -1;
              offset = offset + 2;
            end
          else
            for k = 2:nsectors
              Ai(offset+1, obj.vars.cos_sector.i(k,j-1)) = 1;
              Ai(offset+1, obj.vars.cos_sector.i(k-1:k,j)) = -1;
              Ai(offset+2, obj.vars.sin_sector.i(k,j-1)) = 1;
              Ai(offset+2, obj.vars.sin_sector.i(k-1:k,j)) = -1;
              offset = offset + 2;
            end
          end
        end
        assert(offset == expected_offset);
        obj = obj.addLinearConstraints(Ai, bi, [], []);
      end

      function obj = addInitialSinCosConstraints(obj)
        % Constrain the values of sin and cos for the current poses of the feet
        obj.vars.cos_yaw.lb(1) = cos(obj.seed_plan.footsteps(1).pos(6));
        obj.vars.cos_yaw.ub(1) = cos(obj.seed_plan.footsteps(1).pos(6));
        obj.vars.cos_yaw.lb(2) = cos(obj.seed_plan.footsteps(2).pos(6));
        obj.vars.cos_yaw.ub(2) = cos(obj.seed_plan.footsteps(2).pos(6));
        obj.vars.sin_yaw.lb(3) = sin(obj.seed_plan.footsteps(3).pos(6));
        obj.vars.sin_yaw.ub(3) = sin(obj.seed_plan.footsteps(3).pos(6));
        obj.vars.sin_yaw.lb(4) = sin(obj.seed_plan.footsteps(4).pos(6));
        obj.vars.sin_yaw.ub(4) = sin(obj.seed_plan.footsteps(4).pos(6));
      end

      function plan = getFootstepPlan(obj)
        % Solve the problem if needed and retrieve a footstep plan with the corresponding solution.
        if ~isfield(obj.vars.footsteps, 'value')
          obj = obj.solve();
        end
        steps = zeros(6, obj.nsteps);
        steps(obj.pose_indices, :) = obj.vars.footsteps.value;
        plan = obj.seed_plan;
        for j = 1:obj.nsteps
          plan.footsteps(j).pos = steps(:,j);
        end

        for j = 5:obj.nsteps
          region_ndx = find(obj.vars.region.value(:,j));
          assert(length(region_ndx) == 1, 'Got no (or multiple) region assignments for this footstep. This indicates an infeasibility or bad setup in the mixed-integer program');
          plan.region_order(j) = region_ndx;
        end
        
        if obj.plan_gait
          %obj.vars.desv.value
          %obj.vars.dh.value
          gait_init = int8(obj.vars.seq.value);
          gait_idx  = int8(obj.vars.timming.value);
          idx = max(gait_idx);
          gait_fin = zeros(4,idx);
          for i = 1:4:length(gait_idx)-3
            gait_fin(1,int8(gait_idx(i)))   =  obj.vars.seq.value(i,int8(gait_idx(i)));
            gait_fin(2,int8(gait_idx(i+1))) =  obj.vars.seq.value(i+1,int8(gait_idx(i+1)));
            gait_fin(3,int8(gait_idx(i+2))) =  obj.vars.seq.value(i+2,int8(gait_idx(i+2)));
            gait_fin(4,int8(gait_idx(i+3))) =  obj.vars.seq.value(i+3,int8(gait_idx(i+3)));
          end
          i = 1;
          while 1
            if gait_fin(1:4,i) == zeros(4,1)
              for j = i:length(gait_fin)-1
                gait_fin(1:4,j) = gait_fin(1:4,j+1);
              end
              gait_fin = gait_fin(:,1:end-1);
              i = i - 1;
            end
            if i == size(gait_fin,2)
              break;
            else
              i = i + 1;
            end
          end
          plan.gait = gait_fin;
          plan.timming = gait_idx;
        else
          gait_ex = [1 0 0 0;0 0 1 0;0 0 0 1;0 1 0 0];
          plan.gait = repmat(gait_ex,1,(obj.nsteps-4)/4);
          plan.timming = 1:(obj.nsteps-4);
        end

        %plan = plan.trim_duplicates();
      end
    end
  end

