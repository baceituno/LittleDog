classdef Quad_verifier < Quad_MixedIntegerConvexProgram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
% Trajectory verifier:                                    %
%                                                         %
% Checks if a CoM trajectory satisfies the CWS stability  %
% constraints for a given contact plan.                   %
%                                                         %
% Bernardo Aceituno C. (USB Mechatronics Research Group)  %
% Hongkai Dai          (MIT Robot Locomotion Group / TRI) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  properties(SetAccess = protected)
    nT
    t % a nT x 1 vector
    
    robot_mass;
    gravity;
    hdot % A 3 x nT matrix,
    wrench_disturbance_pt % A 3 x nT matrix, wrench-disturbance_pt(:,i) is the location where the wrench disturbance is applied at time t(i)
    ddr  % A 3 x nT matrix,
    r  % A 3 x nT matrix,
    cws_planner
    parser
  end
  
  methods
    function obj = Quad_verifier(cws_planner, parser)
      obj = obj@Quad_MixedIntegerConvexProgram(false);
      obj.t = cws_planner.t;
      obj.nT = length(cws_planner.t);

      obj.robot_mass = cws_planner.robot_mass;
      obj.gravity = 9.81;
      obj.parser = parser;
      
      obj.wrench_disturbance_pt = cws_planner.wrench_disturbance_pt;
      obj.hdot = cws_planner.vars.dkO.value;
      obj.ddr = cws_planner.vars.ddr.value;
      obj.r = cws_planner.vars.r.value;

      l = length(obj.parser.walking_plan.support_times);
      
      obj = obj.addVariable('F_sum','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('pxF_sum','C',[3,obj.nT],-inf,inf); 
      obj = obj.addVariable('pxF','C',[3,4,obj.nT],-inf,inf); 
      obj = obj.addVariable('F','C',[3,4,obj.nT],-inf,inf);
      obj = obj.addVariable('w','C',[4,4,obj.nT],0,inf);

      display('adding friction cone constraints')
      obj = obj.addFrictionConeconstraint();
      display('adding CoM dynamics constraints')
      obj = obj.addCoMDynamicsconstraint();
      display('adding k0 constraints')
      obj = obj.addk0constraint();
    end
    
    function obj = addFrictionConeconstraint(obj)
      for i = 2:obj.parser.num_supports
        %gets the support interval and time index
        t0 = obj.parser.walking_plan.support_times(i-1);
        t1 = obj.parser.walking_plan.support_times(i);
        t_idx = find(obj.t >= t0 & obj.t < t1);
        %gets the id for the footsteps in contact
        footstep_id = find(obj.parser.footstep_support_interval(1,:) <= i-1 & obj.parser.footstep_support_interval(2,:) > i-1);

        %adds the friction cone constraint for each force
        for j = 1:length(footstep_id)
            %gets cone edges and contact position
            edges_pos = obj.parser.footstep_fc_edges{footstep_id(j)};
            contact_pos = obj.parser.footstep_contact_pos{footstep_id(j)};
            %constraints the force to lie in the friction polyhedral on all timesteps of the support phase
            for k = 1:length(t_idx)
              %constrains in x,y and z
              for ii = 1:3
                Aeq = sparse(1,obj.nv);
                beq = zeros(1,1);
                Aeq(1,obj.vars.F.i(ii,j,t_idx(k))) = 1;
                Aeq(1,obj.vars.w.i(1,j,t_idx(k))) = -edges_pos(ii,1);
                Aeq(1,obj.vars.w.i(2,j,t_idx(k))) = -edges_pos(ii,2);
                Aeq(1,obj.vars.w.i(3,j,t_idx(k))) = -edges_pos(ii,3);
                Aeq(1,obj.vars.w.i(4,j,t_idx(k))) = -edges_pos(ii,4);
                obj = obj.addLinearConstraints([],[],Aeq,beq);
              end
              
              %Gets the cross product between contact position and forces
              Aeq = sparse(3,obj.nv);
              beq = zeros(3,1);

              Aeq(1,obj.vars.pxF.i(1,j,t_idx(k))) = -1;
              Aeq(2,obj.vars.pxF.i(2,j,t_idx(k))) = -1;
              Aeq(3,obj.vars.pxF.i(3,j,t_idx(k))) = -1;
              
              Aeq(1,obj.vars.F.i(3,j,t_idx(k))) = contact_pos(2);
              Aeq(1,obj.vars.F.i(2,j,t_idx(k))) = -contact_pos(3);

              Aeq(2,obj.vars.F.i(1,j,t_idx(k))) = contact_pos(3);
              Aeq(2,obj.vars.F.i(3,j,t_idx(k))) = -contact_pos(1);

              Aeq(3,obj.vars.F.i(2,j,t_idx(k))) = contact_pos(1);
              Aeq(3,obj.vars.F.i(1,j,t_idx(k))) = -contact_pos(2);

              obj = obj.addLinearConstraints([],[],Aeq,beq);       
            end
        end
      end

      for i = 1:obj.nT
        Aeq = sparse(3,obj.nv);
        beq = zeros(3,1);

        %Gets the sum of all forces
        for k = 1:3
          Aeq(k,obj.vars.F_sum.i(k,i)) = 1;
          Aeq(k,obj.vars.F.i(k,1,i)) = -1;
          Aeq(k,obj.vars.F.i(k,2,i)) = -1;
          Aeq(k,obj.vars.F.i(k,3,i)) = -1;
          Aeq(k,obj.vars.F.i(k,4,i)) = -1;
        end

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(3,obj.nv);
        beq = zeros(3,1);

        %Gets the sum of all forces
        for k = 1:3
          Aeq(k,obj.vars.pxF_sum.i(k,i)) = 1;
          Aeq(k,obj.vars.pxF.i(k,1,i)) = -1;
          Aeq(k,obj.vars.pxF.i(k,2,i)) = -1;
          Aeq(k,obj.vars.pxF.i(k,3,i)) = -1;
          Aeq(k,obj.vars.pxF.i(k,4,i)) = -1;
        end

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end

    function obj = addCoMDynamicsconstraint(obj)
      for i=1:obj.nT
        Aeq = sparse(3,obj.nv);
        beq = zeros(3,1);

        Aeq(1,obj.vars.F_sum.i(1,i)) = 1;
        beq(1,1) = obj.robot_mass*obj.ddr(1,i);

        Aeq(2,obj.vars.F_sum.i(2,i)) = 1;
        beq(2,1) = obj.robot_mass*obj.ddr(2,i);

        Aeq(3,obj.vars.F_sum.i(3,i)) = 1;
        beq(3,1) = obj.robot_mass*obj.ddr(3,i) + obj.robot_mass*obj.gravity;

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end

    function obj = addk0constraint(obj)
      for i=1:obj.nT
        Aeq = sparse(3,obj.nv);
        beq = zeros(3,1);
        m = obj.robot_mass;
        Aeq(1,obj.vars.pxF_sum.i(1,i)) = 1;
        beq(1,1) = obj.hdot(1,i) + m*obj.gravity*obj.r(2,i);

        Aeq(2,obj.vars.pxF_sum.i(2,i)) = 1;
        beq(2,1) = obj.hdot(2,i) - m*obj.gravity*obj.r(1,i);

        Aeq(3,obj.vars.pxF_sum.i(3,i)) = 1;
        beq(3,1) = obj.hdot(3,i);

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end
    end

  end
  
end
