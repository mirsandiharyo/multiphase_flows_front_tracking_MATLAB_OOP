% Face class.
classdef Face < handle
    properties(SetAccess = private)
        u
        uOld
        uTemp
        v
        vOld
        vTemp
    end
    
    properties
        forceX
        forceY        
    end
    
    methods
        %% Initialize variables stored at cell Face (liquid is at rest at
        % the beginning).
        function obj = Face(domain)
            % velocity in x-direction
            [obj.u, obj.uOld, obj.uTemp] = deal(zeros(domain.nx+1, domain.ny+2));
            % velocity in y-direction
            [obj.v, obj.vOld, obj.vTemp] = deal(zeros(domain.nx+2, domain.ny+1));
            % forces
            [obj.forceX, obj.forceY] = deal(zeros(domain.nx+2, domain.ny+2));
        end
        
        %% Store old variables for second order scheme.
        function storeOldVariables(obj)
            obj.uOld = obj.u;
            obj.vOld = obj.v;
        end
        
        %% Store second order variables.
        function store2ndOrderVariables(obj)
           obj.u = 0.5*(obj.u+obj.uOld);
           obj.v = 0.5*(obj.v+obj.vOld);
        end
        
        %% Set the forces to zero.
        function initializeForce(obj, domain)
           [obj.forceX, obj.forceY] = deal(zeros(domain.nx+2, domain.ny+2)); 
        end
        
        %% Update the tangential velocity at boundaries.
        function updateWallVelocity(obj, domain)
            % The domain is currently assumed as a box with no-slip boundary
            % condition.
            uSouth = 0;
            uNorth = 0;
            vWest = 0;
            vEast = 0;
            obj.u(1:domain.nx+1,1) = 2*uSouth-obj.u(1:domain.nx+1,2);
            obj.u(1:domain.nx+1,domain.ny+2) = ...
                2*uNorth-obj.u(1:domain.nx+1,domain.ny+1);
            obj.v(1,1:domain.ny+1) = 2*vWest -obj.v(2,1:domain.ny+1);
            obj.v(domain.nx+2,1:domain.ny+1) = ...
                2*vEast -obj.v(domain.nx+1,1:domain.ny+1);
        end
    end
end
