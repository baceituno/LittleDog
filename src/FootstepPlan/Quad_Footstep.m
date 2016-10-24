classdef Quad_Footstep
% A data structure for a single footstep position.
% The pose of the footstep is assumed to be expressed as the [x,y,z,r,p,y] pose
% of the center sole of the foot, which is the location expressed by the Drake
% frame given by frame_id
  properties
    pos
    id
    frame_id
    is_in_contact
    pos_fixed
    terrain_pts
    infeasibility
    walking_params
  end

  methods
    function obj = Quad_Footstep(pos, id, frame_id, is_in_contact, pos_fixed, terrain_pts, infeasibility, walking_params)
      obj.pos = pos;
      obj.id = id;
      obj.frame_id = frame_id;
      obj.is_in_contact = is_in_contact;
      obj.pos_fixed = pos_fixed;
      obj.terrain_pts = terrain_pts;
      obj.infeasibility = infeasibility;
      obj.walking_params = walking_params;
    end
  end
  methods(Static=true)
  end
end




