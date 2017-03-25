classdef Quad_ZMPplanner < Quad_MixedIntegerConvexProgram
  properties(SetAccess = protected)
    nT
    t % a nT x 1 vector
    
    robot_mass;
    gravity;
    
    Ain_cws % A nT x 1 cell
    bin_cws % A nT x 1 cell
    wrench_disturbance_pt % A 3 x nT matrix, wrench-disturbance_pt(:,i) is the location where the wrench disturbance is applied at time t(i)
  end
  
  properties(Access = protected)
    kG_weight_permute
  end
  
  methods
    function obj = Quad_ZMPplanner(robot_mass,parser)
      obj = obj@Quad_MixedIntegerConvexProgram(false);
      
      obj.t = parser.t;
      obj.nT = length(obj.t);
      
      obj.gravity = 9.81;
      
      obj.Ain_cws = parser.Ain_cws;
      obj.bin_cws = parser.bin_cws;

      obj.wrench_disturbance_pt = parser.wrench_disturbance_pt;
      obj = obj.setupVariable();
      obj = obj.integrationConstraint();
    end
  end
  
  methods(Access = protected)
    function obj = setupVariable(obj)
      obj = obj.addVariable('r','C',[3,obj.nT],-inf,inf); 
      obj = obj.addVariable('dr','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('ddr','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('zmp','C',[2,obj.nT],-inf,inf);
      obj = obj.addVariable('dzmp','C',[2,obj.nT],-inf,inf);
      obj = obj.addVariable('com_ellipsoid_slack3','C',[1,obj.nT],1,1);
      obj = obj.addVariable('com_ellipsoid_slack4','C',[3,obj.nT],-inf,inf);
    end
    
    function obj = integrationConstraint(obj)
      dt = diff(obj.t);
      A1_row = reshape(bsxfun(@times,(1:3*(obj.nT-1))',ones(1,4)),[],1);
      A1_col = [(1:3*(obj.nT-1))';3+(1:3*(obj.nT-1))';3*obj.nT+(1:3*(obj.nT-1))';3*obj.nT+3+(1:3*(obj.nT-1))'];
      A1_val = [-ones(3*(obj.nT-1),1);ones(3*(obj.nT-1),1);reshape(bsxfun(@times,-0.5*dt',ones(3,1)),[],1);reshape(bsxfun(@times,-0.5*dt',ones(3,1)),[],1)];
      var_idx = [obj.vars.r.i(:);obj.vars.dr.i(:)];
      Aeq1 = sparse(A1_row,var_idx(A1_col),A1_val,3*(obj.nT-1),obj.nv);
      obj = obj.addLinearConstraints([],[],Aeq1,zeros(3*(obj.nT-1),1));
      
      A2_row = reshape(bsxfun(@times,(1:3*(obj.nT-1))',ones(1,3)),[],1);
      A2_col = [(1:3*(obj.nT-1))';3+(1:3*(obj.nT-1))';3*obj.nT+3+(1:3*(obj.nT-1))'];
      A2_val = [-ones(3*(obj.nT-1),1);ones(3*(obj.nT-1),1);reshape(bsxfun(@times,-dt',ones(3,1)),[],1)];
      var_idx = [obj.vars.dr.i(:);obj.vars.ddr.i(:)];
      Aeq2 = sparse(A2_row,var_idx(A2_col),A2_val,3*(obj.nT-1),obj.nv);
      obj = obj.addLinearConstraints([],[],Aeq2,zeros(3*(obj.nT-1),1));
    end
    
    function obj = addZMPconstraint(obj,com_guess)
      % computes the ZMP at each timestep
      for j = 1:obj.nT
        Aeq = sparse(2,obj.nv);
        beq = zeros(2,1);
        
        Aeq(1,obj.vars.zmp.i(1,j)) = 1;
        Aeq(1,obj.vars.r.i(1,j)) = -1;
        Aeq(1,obj.vars.ddr.i(1,j)) = com_guess/obj.gravity;

        Aeq(2,obj.vars.zmp.i(2,j)) = 1;
        Aeq(2,obj.vars.r.i(2,j)) = -1;
        Aeq(2,obj.vars.ddr.i(2,j)) = com_guess/obj.gravity;

        obj = obj.addLinearConstraints([],[],Aeq,beq);

        Aeq = sparse(2,obj.nv);
        beq = zeros(2,1);
        
        Aeq(1,obj.vars.zmp.i(1,j)) = 1;
        Aeq(1,obj.vars.dzmp.i(1,j)) = -1;
        beq(1,1) = obj.wrench_disturbance_pt(1,j);

        Aeq(2,obj.vars.zmp.i(2,j)) = 1;
        Aeq(2,obj.vars.dzmp.i(2,j)) = -1;
        beq(2,1) = obj.wrench_disturbance_pt(2,j);

        obj = obj.addLinearConstraints([],[],Aeq,beq);
      end

      obj = obj.addVariable('dzmp_norm','C',[1,1],0,inf);
      obj = obj.addConesByIndex([obj.vars.dzmp_norm.i;obj.vars.dzmp.i(:)]);

      c = zeros(obj.nv,1);
      c(obj.vars.dzmp_norm.i) = 1;
      obj = obj.addCost([],c,[]);

      obj = obj.addVariable('accl_norm','C',[1,1],0,inf);
      obj = obj.addConesByIndex([obj.vars.accl_norm.i;obj.vars.ddr.i(:)]);

      c = zeros(obj.nv,1);
      c(obj.vars.accl_norm.i) = 1;
      obj = obj.addCost([],c,[]);
    end
    
    function obj = addCoMEllipsoidBound(obj,A_com,com_guess)
      % add the constraint that |A_com\(com-com_guess)| <= 1
      Aeq3 = sparse(3,obj.nv);
      A_com_inv = inv(A_com);
      for i = 1:obj.nT
        Aeq3((i-1)*3+(1:3),obj.vars.r.i(:,i)) = A_com_inv;
        Aeq3((i-1)*3+(1:3),obj.vars.com_ellipsoid_slack4.i(:,i)) = -eye(3);
      end
      beq3 = reshape(A_com\com_guess,[],1);
      obj = obj.addLinearConstraints([],[],Aeq3,beq3);
      obj = obj.addConesByIndex([obj.vars.com_ellipsoid_slack3.i;obj.vars.com_ellipsoid_slack4.i]);
    end
  end
end