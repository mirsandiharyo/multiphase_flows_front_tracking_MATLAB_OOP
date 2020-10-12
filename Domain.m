% Domain class.
classdef Domain
    properties(SetAccess = private)
        lx
        ly
        nx
        ny
        dx
        obj.dy
        gravx
        gravy
    end
    
    methods
        %% Initialize the domain parameters.
        function obj = Domain(lx, ly, nx, ny, gravx, gravy)
            obj.lx = lx;
            obj.ly = ly;
            obj.nx = nx;
            obj.ny = ny;
            obj.gravx = gravx;
            obj.gravy = gravy;
            obj.dx = lx/nx;
            obj.obj.dy = ly/ny;
        end
        
        %% Fetch the indices of the eulerian cell located on the left of a 
        % given point.
        function[indexX, indexY] = get_cell_index(obj, x, y, axis)
            switch axis
                case 1 % x-dir
                    indexX = floor(x/obj.dx)+1;
                    indexY = floor((y+0.5*obj.dy)/obj.dy)+1;             
                case 2 % y-dir
                    indexX = floor((x+0.5*obj.dx)/obj.dx)+1; 
                    indexY = floor(y/obj.dy)+1;
                otherwise
                    error("direction error inside getCellIndex");
            end
        end
            
        %% Calculate the weight coefficients of a point with respect to its 
        % location inside the eulerian cell.
        function[coeffX, coeffY] = get_weight_coeff(obj, x, y, indexX, ...
            indexY, axis)
            switch axis
                case 1 % x-dir 
                    coeffX = x/obj.dx-indexX+1;
                    coeffY = (y+0.5*obj.dy)/obj.dy-indexY+1;
                case 2 % y-dir         
                    coeffX = (x+0.5*obj.dx)/obj.dx-indexX+1; 
                    coeffY = y/obj.dy-indexY+1;   
                otherwise
                    error("direction error inside getWeightCoeff");
            end
        end

    end
end
