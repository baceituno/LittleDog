classdef Quadruped < LeggedRobot
% Interface class for an Quadruped robot. Being a quadruped currently affords
% footstep planning and ZMP-based walking planning, but these capabilities
% will expand in the future.
% In order to function as an Quadruped, a robot's URDF must be tagged with six
% Drake-style <frame> tags indicating the centers of the soles of its feet.

  properties
    foot_body_id
    foot_frame_id
    inner_foot_shape
    default_body_collision_slices % see getBodyCollisionSlices() below
    body_collision_slice_heights % see getBodyCollisionSlices() below

  end

  properties(Abstract)
    default_footstep_params
    default_walking_params

    base_link;

    Foot1;
    Foot2;
    Foot3;
    Foot4;
  end

  methods
    function obj = Quadruped(Foot1, Foot2, Foot3, Foot4)
      % Construct a Quadruped by identifying the Drake frame's corresponding
      % to the soles of its feet.
      if nargin < 4
        % use littledog defaults
        Foot1 = 'front_left_foot_center';
        Foot2 = 'front_right_foot_center';
        Foot3 = 'back_left_foot_center';
        Foot4 = 'back_right_foot_center';
      end

      obj = obj@LeggedRobot();

      Foot1_frame_id = obj.parseBodyOrFrameID(Foot1);
      Foot1_body_id = obj.getFrame(Foot1_frame_id).body_ind;
      Foot2_frame_id = obj.parseBodyOrFrameID(Foot2);
      Foot2_body_id = obj.getFrame(Foot2_frame_id).body_ind;
      Foot3_frame_id = obj.parseBodyOrFrameID(Foot3);
      Foot3_body_id = obj.getFrame(Foot3_frame_id).body_ind;
      Foot4_frame_id = obj.parseBodyOrFrameID(Foot4);
      Foot4_body_id = obj.getFrame(Foot4_frame_id).body_ind;
      obj.foot_body_id = struct('Foot1', Foot1_body_id, 'Foot2', Foot2_body_id, 'Foot3', Foot3_body_id,...
                                'Foot4', Foot4_body_id);
      
      obj.foot_frame_id = struct('Foot1', Foot1_frame_id, 'Foot2', Foot2_frame_id, 'Foot3', Foot3_frame_id,...
                                 'Foot4', Foot4_frame_id);
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
      
      weights = struct('relative', [1;1;1;0;0;0],...
                       'relative_final', [10;10;10;0;0;0],...
                       'goal', [1000;1000;0;0;0;100]);
    end

    function params = applyDefaultFootstepParams(obj, params)
      params = applyDefaults(params, obj.default_footstep_params);
    end

    function foot_center = feetPosition(obj, q0)
      % Convenient way to find the poses of the center soles of the feet given a
      % configuration vector q0

      sizecheck(q0,[obj.getNumPositions,1]);

      kinsol = doKinematics(obj,q0);

      Foot1 = forwardKin(obj,kinsol,obj.foot_frame_id.Foot1,[0;0;0], 1);
      Foot2 = forwardKin(obj,kinsol,obj.foot_frame_id.Foot2,[0;0;0], 1);
      Foot3 = forwardKin(obj,kinsol,obj.foot_frame_id.Foot3,[0;0;0], 1);
      Foot4 = forwardKin(obj,kinsol,obj.foot_frame_id.Foot4,[0;0;0], 1);
      
      Foot1(6,1) = q0(6,1);
      Foot2(6,1) = q0(6,1);
      Foot3(6,1) = q0(6,1);
      Foot4(6,1) = q0(6,1);

      foot_center = struct('Foot1', Foot1, 'Foot2', Foot2, 'Foot3', Foot3, 'Foot4', Foot4);
    end

    function foot_min_z = getFootHeight(obj,q)
      % Get the height in world coordinates of the lower of the robot's foot soles
      foot_center = obj.feetPosition(q);
      foot_min_z = min(foot_center.Foot1(3), foot_center.Foot2(3), foot_center.Foot3(3),...
                       foot_center.Foot4(3));
    end
  
    function fc = getFootContacts(obj, q)
      % For a given configuration of the Quadruped, determine whether each foot is
      % in contact with the terrain. 
      % @param q a robot configuration vector
      % @retval fc a logical vector of length 6. If fc(i) is true, then foot i is in contact
      [phiC,~,~,~,idxA,idxB] = obj.collisionDetect(q,false);
      within_thresh = phiC < 0.002;
      contact_pairs = [idxA(within_thresh); idxB(within_thresh)];
      
      % The following would be faster but would require us to have
      % heightmaps in Bullet
      %[~,~,idxA,idxB] = obj.r_control.allCollisions(x(1:obj.nq_control));
      %contact_pairs = [idxA; idxB];
      foot_indices = [obj.foot_body_id.Foot1, obj.foot_body_id.Foot2, obj.foot_body_id.Foot3,...
                      obj.foot_body_id.Foot4];
      fc = any(bsxfun(@eq, contact_pairs(:), foot_indices),1)';
    end

    function collision_model = getFootstepPlanningCollisionModel(obj, varargin)
      % Get a simple collision model for the Atlas robot to enable
      % (limited) whole-body awareness during footstep planning. The
      % collision model is represented as a set of points which define
      % the shape of the foot and a set of slices which define the 
      % bounding boxes of the legs, torso, and arms. The collision model
      % of the foot used here may actually be smaller than the robot's
      % feet in order to allow the toes or heels to hang over edges. 
      % @param q (optional) if provided, use the given robot configuration 
      %          to compute the collision volumes. Otherwise use hard-coded
      %          values derived from a typical walking posture.
      % @retval collision_model an IRIS CollisionModel object with fields
      %         foot and body. The body field has subfields z and xy, where
      %         xy is of shape [3, N, length(z)]. Each page of xy(:,:,j)
      %         represents the bounds of the robot from z(j) to z(j+1) (or inf).
      checkDependency('iris');
      foot_shape = obj.getInnerFootShape();
      slices = obj.getBodyCollisionSlices(varargin{:});
      collision_model = iris.terrain_grid.CollisionModel(foot_shape, slices);
    end

    function shape = getInnerFootShape(obj)
      % Get an inner approximation of the foot shape as a set of points in xy. The convention here is that this entire shape must be supported by the terrain. So, making this shape smaller than the actual foot will allow the edge of the foot to hang over edges when the footstep planner is run. 
      % By convention, x points forward (from heel to toe), y points left, and the foot shape is symmetric about y=0.
      if ~isempty(obj.inner_foot_shape)
        shape = obj.inner_foot_shape;
      else
        shape = [-0.01, -0.01, 0.01,  0.01;
                  0.01, -0.01, 0.01, -0.01];
      end
    end
  end
end
