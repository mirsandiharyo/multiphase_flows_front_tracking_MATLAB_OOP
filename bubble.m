% Bubble class.
classdef bubble
    properties(SetAccess = private)
       center_x
       center_y
       radius
       point
       x
       y
       x_old
       y_old
    end
    
    methods
        % Initialize the bubble.
        function obj = bubble(center_x, center_y, radius, point)
            obj.center_x = center_x;
            obj.center_y = center_y;
            obj.radius = radius;
            obj.point = point;
            [obj.x, obj.y, obj.x_old, obj.y_old] = ...
            deal(zeros(1,obj.point+2));
        end
    end
end
