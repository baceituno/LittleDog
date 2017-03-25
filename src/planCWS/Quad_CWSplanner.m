classdef Quad_CWSplanner < Quad_MixedIntegerConvexProgram
  properties(SetAccess = protected)
    nT
    t % a nT x 1 vector
    
    robot_mass;
    gravity;
    
    Ain_cws % A nT x 1 cell
    bin_cws % A nT x 1 cell
    wrench_disturbance_pt % A 3 x nT matrix, wrench-disturbance_pt(:,i) is the location where the wrench disturbance is applied at time t(i)
    
    Q_w
    
    centroidal_angular_momentum_norm_weight % A 3 x 1 vector. The norm of the centroidal angular momentum kG is defined as weight'*abs(kG)
  end
  
  properties(Access = protected)
    kG_weight_permute
  end
  
  methods
    function obj = Quad_CWSplanner(robot_mass,parser,Q_w,centroidal_angular_momentum_norm_weight)
      obj = obj@Quad_MixedIntegerConvexProgram(false);
      
      obj.t = parser.t;
      obj.nT = length(obj.t);
      if(numel(robot_mass) ~= 1 || robot_mass <= 0)
          error('robot_mass should be a positive scalar');
      end
      obj.robot_mass = robot_mass;
      obj.gravity = 9.81;
      
      if(~iscell(parser.Ain_cws) || ~iscell(parser.bin_cws) || length(parser.Ain_cws)~=obj.nT || length(parser.bin_cws) ~= obj.nT)
        error('Ain_cws and bin_cws do not have the right size');
      end
      obj.Ain_cws = parser.Ain_cws;
      obj.bin_cws = parser.bin_cws;
      
      if(any(size(centroidal_angular_momentum_norm_weight) ~= [3,1]) || any(centroidal_angular_momentum_norm_weight < 0))
          error('centroidal_angular_momentum_norm_weight should be a 3 x 1 non-negative vector');
      end
      obj.centroidal_angular_momentum_norm_weight = centroidal_angular_momentum_norm_weight;
      obj.kG_weight_permute = repmat(obj.centroidal_angular_momentum_norm_weight,1,8).*[1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 1 1 -1 -1; 1 -1 1 -1 1 -1 1 -1];
      
      if(any(size(parser.wrench_disturbance_pt) ~= [3,obj.nT]))
          error('wrench_disturbance_pt should be of size 3 x %d',obj.nT);
      end
      obj.wrench_disturbance_pt = parser.wrench_disturbance_pt;
      if(any(size(Q_w) ~= [6,6]))
        error('Qw should be of size 6 x 6');
      end
      Q_w = (Q_w+Q_w')/2;
      if(any(eig(Q_w) < 0))
        error('Qw should be psd');
      end
      obj.Q_w = Q_w;
      obj = obj.setupVariable();
      obj = obj.integrationConstraint();
      obj = obj.CWSconstraint();
    end
    
    function obj = addCWSmarginLowerBound(obj,cws_margin_lb)
      obj.vars.cws_margin.lb = max([obj.vars.cws_margin.lb;cws_margin_lb*ones(1,obj.nT)],[],1);
    end
    
    function obj = addCentroidalAngularMomentumCost(obj)
      c = zeros(obj.nv,1);
      c(obj.vars.kG_norm_ub.i) = 1;
      obj = obj.addCost([],c,[]);
    end
    
    function obj = addCOMaccelCost(obj,Q)
      sizecheck(Q,[3,3]);
      Q_row = reshape(repmat(obj.vars.ddr.i,3,1),[],1);
      Q_col = reshape(repmat(reshape(obj.vars.ddr.i,1,[]),3,1),[],1);
      Q_val = reshape(repmat(Q(:),1,obj.nT),[],1);
      Q_all = sparse(Q_row,Q_col,Q_val,obj.nv,obj.nv);
      obj = obj.addCost(Q_all,[],[]);
    end
    
    function obj = addCWSmarginCost(obj,weight)
      if(numel(weight) ~= 1 || weight<=0)
        error('weight should be a positive scalar');
      end
      c = zeros(obj.nv,1);
      c(obj.vars.cws_margin.i(2:end)) = -weight;
      obj = obj.addCost([],c,[]);
    end
  end
  
  methods(Access = protected)
    function obj = setupVariable(obj)
      obj = obj.addVariable('r','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('dr','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('ddr','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('kO','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('dkO','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('cws_margin','C',[1,obj.nT],0,inf);
      obj = obj.addVariable('kG_norm_ub','C',[1,obj.nT],0,inf);
      obj = obj.addVariable('com_ellipsoid_slack1','C',[8,obj.nT],0,inf);
      obj = obj.addVariable('com_ellipsoid_slack2','C',[24,obj.nT],-inf,inf);
      obj = obj.addVariable('com_ellipsoid_slack3','C',[1,obj.nT],1,1);
      obj = obj.addVariable('com_ellipsoid_slack4','C',[3,obj.nT],-inf,inf);
      obj = obj.addVariable('com_accel_norm','C',[1,1],0,inf);
      obj = obj.addConesByIndex([obj.vars.com_accel_norm.i;obj.vars.ddr.i(:)]);
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
      
      var_idx = [obj.vars.kO.i(:);obj.vars.dkO.i(:)];
      Aeq3 = sparse(A2_row,var_idx(A2_col),A2_val,3*(obj.nT-1),obj.nv);
      obj = obj.addLinearConstraints([],[],Aeq3,zeros(3*(obj.nT-1),1));
    end
    
    function obj = CWSconstraint(obj)
      % compute the CWS margin
      iA = cell(obj.nT,1);
      jA = cell(obj.nT,1);
      Aval = cell(obj.nT,1);
      ub = cell(obj.nT,1);
      num_A_rows = 0;
      x_ind = cell(obj.nT,1);
      for i = 1:obj.nT
        Ai = sparse(length(obj.bin_cws{i}),10);

        ub{i} = obj.bin_cws{i} + obj.Ain_cws{i}(:,1:3)*[0;0;-obj.robot_mass*obj.gravity];
        Ai(:,1:6) = [obj.Ain_cws{i}(:,1:3)*obj.robot_mass obj.Ain_cws{i}(:,4:6)];
        Ai(:,7:9) = -obj.Ain_cws{i}(:,4:6)*crossSkewSymMat([0;0;obj.robot_mass*obj.gravity]);
        T = [eye(3) zeros(3);crossSkewSymMat(obj.wrench_disturbance_pt(:,i)) eye(3)];
        Ai(:,10) = sqrt(sum((obj.Ain_cws{i}*T/obj.Q_w).*(obj.Ain_cws{i}*T),2));
        
        [iA{i},jA{i}] = find(Ai);
        Aval{i} = Ai(sub2ind(size(Ai),iA{i},jA{i}));
        iA{i} = iA{i} + num_A_rows;
        jA{i} = jA{i} + 10*(i-1);
        
        x_ind{i} = [obj.vars.ddr.i(:,i);obj.vars.dkO.i(:,i);obj.vars.r.i(:,i);obj.vars.cws_margin.i(i)];
        num_A_rows = num_A_rows+size(Ai,1);
      end
      x_ind_all = vertcat(x_ind{:});
      Ain = sparse(vertcat(iA{:}),x_ind_all(vertcat(jA{:})),vertcat(Aval{:}),num_A_rows,obj.nv);
      bin = vertcat(ub{:});
      obj = obj.addLinearConstraints(Ain,bin,[],[]);
    end
    
    function obj = addCentroidalAngularMomentumBndCOMellipsoid(obj,kG_norm_ub_idx,com_guess,A_com)
      % add the constraint that
      % kG_norm_ub-weight*(kO+m*cross(dr,com_guess))>=m*|A_com'*cross(weight,com_vel)|
      % first add com_ellipsoid_slack1 =
      % kG_norm_ub-weight'*(kO+m*cross(dr,com_guess)
      iA1 = cell(obj.nT,1);
      jA1 = cell(obj.nT,1);
      Aval1 = cell(obj.nT,1);
      num_A_rows1 = 0;
      iA2 = cell(obj.nT,1);
      jA2 = cell(obj.nT,1);
      Aval2 = cell(obj.nT,1);
      num_A_rows2 = 0;
      for i = 1:obj.nT
        A1 = sparse(8,obj.nv);
        A1(:,kG_norm_ub_idx(i)) = ones(8,1);
        A1(:,obj.vars.kO.i(:,i)) = -obj.kG_weight_permute';
        A1(:,obj.vars.dr.i(:,i)) = obj.kG_weight_permute'*obj.robot_mass*crossSkewSymMat(com_guess(:,i));
        A1(:,obj.vars.com_ellipsoid_slack1.i(:,i)) = -eye(8);
        [iA1{i},jA1{i},Aval1{i}] = find(A1);
        iA1{i} = iA1{i}+num_A_rows1;
        num_A_rows1 = num_A_rows1 + 8;
        
        % now add the com_ellipsoid_slack2 =
        % m*A_com'*crossSkewSymMat(weight)*dr
        A2 = sparse(24,27);
        for j = 1:8
          A2((j-1)*3+(1:3),(j-1)*3+(1:3))= eye(3);
          A2((j-1)*3+(1:3),24+(1:3)) = -obj.robot_mass*A_com'*crossSkewSymMat(obj.kG_weight_permute(:,j));          
        end
        [iA2{i},jA2{i},Aval2{i}] = find(A2);
        iA2{i} = iA2{i}+num_A_rows2;
        var_idx = [obj.vars.com_ellipsoid_slack2.i(:,i);obj.vars.dr.i(:,i)];
        jA2{i} = var_idx(jA2{i});
        num_A_rows2 = num_A_rows2 + 24;
      end
      Aeq1 = sparse(vertcat(iA1{:}),vertcat(jA1{:}),vertcat(Aval1{:}),num_A_rows1,obj.nv);
      beq1 = zeros(num_A_rows1,1);
      obj = obj.addLinearConstraints([],[],Aeq1,beq1);
      Aeq2 = sparse(vertcat(iA2{:}),vertcat(jA2{:}),vertcat(Aval2{:}),num_A_rows2,obj.nv);
      beq2 = zeros(num_A_rows2,1);
      obj = obj.addLinearConstraints([],[],Aeq2,beq2);
      obj = obj.addConesByIndex([reshape(obj.vars.com_ellipsoid_slack1.i,1,[]);reshape(obj.vars.com_ellipsoid_slack2.i,3,[])]);
      
      % add the constraint that |A_com\(com-com_guess)|<=1
      Aeq3 = sparse(3,obj.nv);
      A_com_inv = inv(A_com);
      obj.nT
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