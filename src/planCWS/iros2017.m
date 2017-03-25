classdef iros2017 < Quad_CWSplanner
  properties(SetAccess = protected)
    com_guess
    A_com
  end
  
  methods
    function obj = iros2017(robot_mass,parser,Q_w,A_com,centroidal_angular_momentum_norm_weight)
      obj = obj@Quad_CWSplanner(robot_mass,parser,Q_w,centroidal_angular_momentum_norm_weight);
      if(any(size(parser.com_guess)~=[3,obj.nT]))
        error('com_guess should be a 3 x %d matrix',obj.nT);
      end
      obj.com_guess = parser.com_guess;
      if(any(size(A_com)~=[3,3]))
        error('A_com should be a 3 x 3 matrix');
      end
      obj.A_com = A_com;
      display('adding kG L1 norm bound constraint');
      obj = obj.addCentroidalAngularMomentumBndCOMellipsoid(obj.vars.kG_norm_ub.i,obj.com_guess,A_com);
      obj = obj.addCWSmarginCost(0.01);
      obj = obj.addCentroidalAngularMomentumCost();

      % adds cost to the com acceleration
      c = zeros(obj.nv,1);
      c(obj.vars.com_accel_norm.i) = 1e-4;
      obj = obj.addCost([],c,[]);
    end
  end
end