% Center class.
classdef Center
    properties(SetAccess = private)
        x
        y
        pres
    end
    
    methods
        % Initialize variables stored at cell center.
        function obj = Center(domain)
            % set the grid
            obj.x = linspace(-0.5, domain.nx+2-1.5, domain.nx+2)*domain.dx;
            obj.y = linspace(-0.5, domain.ny+2-1.5, domain.ny+2)*domain.dy;
            % pressure
            obj.pres = zeros(domain.nx+2, domain.ny+2);
        end
    end
end
