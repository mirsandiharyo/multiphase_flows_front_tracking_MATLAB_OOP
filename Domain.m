% Domain class.
classdef Domain
    properties(SetAccess = private)
        lx
        ly
        nx
        ny
        gravx
        gravy
    end
    
    methods
        % Initialize the domain parameters.
        function obj = Domain(lx, ly, nx, ny, gravx, gravy)
            obj.lx = lx;
            obj.ly = ly;
            obj.nx = nx;
            obj.ny = ny;
            obj.gravx = gravx;
            obj.gravy = gravy;
        end
    end
end
