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
        %% Initialize the bubble.
        function obj = Bubble(centerX, centerY, radius, point)
            obj.centerX = centerX;
            obj.centerY = centerY;
            obj.radius = radius;
            obj.point = point;
            [obj.x, obj.y, obj.xOld, obj.yOld] = ...
            deal(zeros(1,obj.point+2));
        end
        
        %% Determine the location of the initial spherical bubble.
        function initializeFront(obj)
            for i=1:obj.point+2
                obj.x(i) = obj.centerX-obj.radius*sin(2.0*pi*(i-1)/(obj.point));
                obj.y(i) = obj.centerY+obj.radius*cos(2.0*pi*(i-1)/(obj.point));
            end
        end
        
        %% Store old variables for second order scheme.
        function storeOldVariables(obj)
            obj.xOld = obj.x;
            obj.yOld = obj.y;
        end
        
        %% Store second order variables.
        function store2ndOrderVariables(obj)
            obj.x = 0.5*(obj.x+obj.xOld);
            obj.y = 0.5*(obj.y+obj.yOld);            
        end
        
        %% Restructure the front to maintain the quality of the interface.
        function restructure_front(obj, domain)
            obj.xOld = obj.x;
            obj.yOld = obj.y;
            j = 1;
            for i=2:obj.point+1
                % check the distance
                dst = sqrt(((obj.xOld(i)-obj.x(j))/domain.dx)^2 + ...
                           ((obj.yOld(i)-obj.y(j))/domain.dy)^2);        
                if (dst > 0.5) % too big
                    % add marker points
                    j = j+1;
                    obj.x(j) = 0.5*(obj.xOld(i)+obj.x(j-1));
                    obj.y(j) = 0.5*(obj.yOld(i)+obj.y(j-1));
                    j = j+1;
                    obj.x(j) = obj.xOld(i);
                    obj.y(j) = obj.yOld(i);
                elseif (dst < 0.25) % too small
                    % do nothing  
                else
                    j = j+1;
                    obj.x(j) = obj.xOld(i);
                    obj.y(j) = obj.yOld(i);
                end
            end
            obj.point = j-1;
            obj.x(1) = obj.x(obj.point+1);
            obj.y(1) = obj.y(obj.point+1);
            obj.x(obj.point+2) = obj.x(2);
            obj.y(obj.point+2) = obj.y(2);
        end
        
    end
    
    methods (Static)
        %% Distribute a value from a lagrangian point to neighboring eulerian cells.
        function cell = distributeLagrangianToEulerian(domain, cell, x, y, ...
                value, axis)
            % assign the grid size
            switch axis
                case 1 % x-dir
                    d1 = domain.dx;
                    d2 = domain.dy;    
                case 2 % y-dir
                    d1 = domain.dy;
                    d2 = domain.dx;   
                otherwise
                    error("direction error inside distributeLagrangianToEulerian");
            end
            % get the eulerian cell indices
            [indexX, indexY] = getCellIndex(x, y, domain.dx, domain.dy, axis); 
            % calculate the weighing coefficients
            [coeffX, coeffY] = getWeightCoeff(x, y, domain.dx, domain.dy, ...
                indexX, indexY, axis); 

            % distribute the force to the surrounding eulerian cells
            cell(indexX  ,indexY  ) = cell(indexX  ,indexY  ) + ...
                (1.0-coeffX)*(1.0-coeffY)*value/d1/d2;
            cell(indexX+1,indexY  ) = cell(indexX+1,indexY  ) + ...
                coeffX*(1.0-coeffY)*value/d1/d2;
            cell(indexX  ,indexY+1) = cell(indexX  ,indexY+1) + ... 
                (1.0-coeffX)*coeffY*value/d1/d2;      
            cell(indexX+1,indexY+1) = cell(indexX+1,indexY+1) + ...
                coeffX*coeffY*value/d1/d2;
        end
    end
end
