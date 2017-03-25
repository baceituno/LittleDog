classdef testIROS < Quad_ZMPplanner
  properties(SetAccess = protected)
    com_guess
    A_com
    Q_w
  end
  
  methods
    function obj = testIROS(robot_mass,parser,A_com, com_height_guess)
      obj = obj@Quad_ZMPplanner(robot_mass,parser);
      if(any(size(parser.com_guess)~=[3,obj.nT]))
        error('com_guess should be a 3 x %d matrix',obj.nT);
      end
      obj.Q_w  = parser.Q_w;
      obj.com_guess = parser.com_guess;
      if(any(size(A_com)~=[3,3]))
        error('A_com should be a 3 x 3 matrix');
      end
      obj.A_com = A_com;
      display('adding kG L1 norm bound constraint');
      obj = obj.addZMPconstraint(com_height_guess);
      obj = obj.addCoMEllipsoidBound(A_com,obj.com_guess);
    end

    function cws_margin = CWS(obj, com_traj)
        for i = 1:obj.nT
          Ain = obj.Ain_cws{i};
          bin = obj.bin_cws{i};

          T_w = [eye(3) zeros(3);crossSkewSymMat(obj.wrench_disturbance_pt(:,i)) eye(3)];
          prog = spotsosprog();
          [prog,central_ray_i] = prog.newFree(6,1);
          prog = prog.withLor([1;central_ray_i]);
          [prog,dist_i] = prog.newPos(1,1);
          prog = prog.withPos(bin-Ain*central_ray_i-sqrt(sum((Ain*T_w/obj.Q_w).*(Ain*T_w),2))*dist_i);
          sol = prog.minimize(-dist_i,@spot_mosek);
          dist = double(sol.eval(dist_i));
          
          cws_margin(i) = dist;
        end
    end
  end
end
