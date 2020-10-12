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

        %% Advect the location of marker points using the interpolated velocity field.
        function updateFrontLocation(obj, face, param, domain)
            % interpolate the velocity from the eulerian grid to the 
            % location of marker point   
            [uX, uY] = deal(zeros(1,obj.point+2));
            for i=2:obj.point+1
                % interpolate velocity in x-direction
                uX(i) = obj.interpolateVelocity(domain, face.u, obj.x(i), ...
                    obj.y(i), 1);
                % interpolate velocity in y-direction
                uY(i) = obj.interpolateVelocity(domain, face.v, obj.x(i), ...
                    obj.y(i), 2);
            end

            % advect the marker point 
            for i=2:obj.point+1
                obj.x(i) = obj.x(i)+param.dt*uX(i);
                obj.y(i) = obj.y(i)+param.dt*uY(i);
            end
            obj.x(1) = obj.x(obj.point+1);
            obj.y(1) = obj.y(obj.point+1);
            obj.x(obj.point+2) = obj.x(2);
            obj.y(obj.point+2) = obj.y(2);              
        end
    
        %% Calculate the surface tension force on the lagrangian grid and 
        % distribute it to the surrounding eulerian grid cells.
        function calculateSurfaceTension(obj, domain, fluidProp, face)
            % initialize the variables to store the tangent vector
            [tanX, tanY] = deal(zeros(obj.point+2, obj.point+2));

            % calculate the tangent vector
            for i=1:obj.point+1
                dist = sqrt((obj.x(i+1)-obj.x(i))^2 + ...
                            (obj.y(i+1)-obj.y(i))^2);
                tanX(i) = (obj.x(i+1)-obj.x(i))/dist;
                tanY(i) = (obj.y(i+1)-obj.y(i))/dist;
            end
            tanX(obj.point+2) = tanX(2);
            tanY(obj.point+2) = tanY(2);

            % distribute the surface tension force to the eulerian grid
            for i=2:obj.point+1
                % force in x-direction
                forceX = fluidProp.sigma*(tanX(i)-tanX(i-1));
                face.forceX = obj.distributeLagrangianToEulerian(domain, ...
                    face.forceX, obj.x(i), obj.y(i), forceX, 1);

                % force in y-direction
                forceY = fluidProp.sigma*(tanY(i)-tanY(i-1));
                face.forceY = obj.distributeLagrangianToEulerian(domain, ...
                    face.forceY, obj.x(i), obj.y(i), forceY, 2);
            end   
        end
    end
    
    methods (Static)
        %% Distribute a value from a lagrangian point to neighboring 
        % eulerian cells.
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
            [indexX, indexY] = domain.getCellIndex(x, y, axis); 
            % calculate the weighing coefficients
            [coeffX, coeffY] = domain.getWeightCoeff(x, y, indexX, indexY, axis); 

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
        
        %% Interpolate velocities located on eulerian cells to a lagrangian 
        % marker point.
        function vel = interpolateVelocity(domain, faceVel, x, y, axis)
            % get the eulerian cell index
            [indexX, indexY] = domain.getCellIndex(x, y, axis); 
            % calculate the weighing coefficient 
            [coeffX, coeffY] = domain.getWeightCoeff(x, y, indexX, indexY, axis); 
            % interpolate the surrounding velocities to the marker location
            vel = (1.0-coeffX)*(1.0-coeffY)*faceVel(indexX  ,indexY  )+ ...
                       coeffX *(1.0-coeffY)*faceVel(indexX+1,indexY  )+ ...
                  (1.0-coeffX)*     coeffY *faceVel(indexX  ,indexY+1)+ ...
                       coeffX *     coeffY *faceVel(indexX+1,indexY+1);            
        end
        
    end
end
