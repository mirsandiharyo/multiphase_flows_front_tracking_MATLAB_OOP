% Face class.
classdef Face
    properties(SetAccess = private)
        u
        uOld
        uTemp
        v
        vOld
        vTemp
        forceX
        forceY
    end
    
    methods
        % Initialize variables stored at cell Face (liquid is at rest at
        % the beginning).
        function obj = Face(domain)
            % velocity in x-direction
            [obj.u, obj.uOld, obj.uTemp] = deal(zeros(domain.nx+1, domain.ny+2));
            % velocity in y-direction
            [obj.v, obj.vOld, obj.vTemp] = deal(zeros(domain.nx+2, domain.ny+1));
            % forces
            [obj.forceX, obj.forceY] = deal(zeros(domain.nx+2, domain.ny+2));
        end
    end
end