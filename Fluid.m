% Fluid class.
classdef Fluid
    properties(SetAccess = private)
        rho
        rhoOld
        mu
        muOld
    end
    
    methods
        % Initialize the density and viscosity fields using the properties
        % from continuous phase.
        function obj = Fluid(domain, fluidProp)
            [obj.rho, obj.rhoOld] = ...
            deal(zeros(domain.nx+2, domain.ny+2)+fluidProp(2).contRho);
            [obj.mu, obj.muOld] =  ...
            deal(zeros(domain.nx+2, domain.ny+2)+fluidProp(2).contMu);
        end
        
        % Set the fluid properties inside the discrete phase with an initial
        % spherical shape. 
        function initializeDomain(domain, center, bubbleList, fluidProp)
            for i=2:domain.nx+1
               for j=2:domain.ny+1
                   for n=1:length(bubbleList)
                      if ((center.x(i)-bubbleList{n}.centerX)^2 + ...
                          (center.y(j)-bubbleList{n}.centerY)^2 < ...
                          bubbleList{n}.radius^2)
                          obj.rho(i,j) = fluidProp(1).dispRho;
                          obj.mu(i,j)  = fluidProp(1).dispMu;
                      end
                   end
               end
            end            
        end
        
    end
end
