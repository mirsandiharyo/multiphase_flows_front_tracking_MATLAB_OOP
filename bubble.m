% Bubble class.
classdef Bubble < handle
    properties(SetAccess = private)
       centerX
       centerY
       radius
       point
       x
       y
       xOld
       yOld
    end
    
    methods
        % Initialize the bubble.
        function obj = Bubble(centerX, centerY, radius, point)
            obj.centerX = centerX;
            obj.centerY = centerY;
            obj.radius = radius;
            obj.point = point;
            [obj.x, obj.y, obj.xOld, obj.yOld] = ...
            deal(zeros(1,obj.point+2));
        end
        
        % Determine the location of the initial spherical bubble.
        function initializeFront(obj)
            for i=1:obj.point+2
                obj.x(i) = obj.centerX-obj.radius*sin(2.0*pi*(i-1)/(obj.point));
                obj.y(i) = obj.centerY+obj.radius*cos(2.0*pi*(i-1)/(obj.point));
            end
        end
    end
end
