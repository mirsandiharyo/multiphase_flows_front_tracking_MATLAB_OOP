% Fluid class.
classdef Fluid < handle
    properties(SetAccess = private)
        rho
        rhoOld
        mu
        muOld
    end
    
    methods
        %% Initialize the density and viscosity fields using the properties
        % from continuous phase.
        function obj = Fluid(domain, fluidProp)
            [obj.rho, obj.rhoOld] = ...
            deal(zeros(domain.nx+2, domain.ny+2)+fluidProp.contRho);
            [obj.mu, obj.muOld] =  ...
            deal(zeros(domain.nx+2, domain.ny+2)+fluidProp.contMu);
        end
        
        %% Set the fluid properties inside the discrete phase with an initial
        % spherical shape. 
        function initializeDomain(obj, domain, center, bubbleList, fluidProp)
            for i=2:domain.nx+1
               for j=2:domain.ny+1
                   for n=1:length(bubbleList)
                      if ((center.x(i)-bubbleList{n}.centerX)^2 + ...
                          (center.y(j)-bubbleList{n}.centerY)^2 < ...
                          bubbleList{n}.radius^2)
                          obj.rho(i,j) = fluidProp.dispRho;
                          obj.mu(i,j)  = fluidProp.dispMu;
                      end
                   end
               end
            end            
        end
        
        %% Store old variables for second order scheme.
        function storeOldVariables(obj)
            obj.rhoOld = obj.rho;
            obj.muOld = obj.mu;
        end
    end
end
