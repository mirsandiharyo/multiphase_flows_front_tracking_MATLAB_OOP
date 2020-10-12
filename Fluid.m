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
        
        %% Store second order variables.
        function store2ndOrderVariables(obj)
           obj.rho = 0.5*(obj.rho+obj.rhoOld);
           obj.mu = 0.5*(obj.mu+obj.muOld);
        end
        
        %% Update the density field using the density jump at the lagrangian 
        % interface. Linear averaging is used to get the value at each cell.
        function updateDensity(obj, param, domain, bubbleList, fluidProp)
            % initialize the variables to store the density jump
            [faceX, faceY] = deal(zeros(domain.nx+2, domain.ny+2));
            
            % distribute the density jump to the eulerian grid
            for n=1:length(bubbleList)
                for i=2:bubbleList{n}.point+1
                    % density jump in x-direction
                    forceX = -0.5*(bubbleList{n}.y(i+1)- ...
                                   bubbleList{n}.y(i-1))* ...
                                  (fluidProp.dispRho-fluidProp.contRho);  
                    faceX = bubbleList{n}.distributeLagrangianToEulerian( ...
                            domain, faceX, bubbleList{n}.x(i), ...
                            bubbleList{n}.y(i), forceX, 1);

                    % density jump in y-direction
                    forceY = 0.5*(bubbleList{n}.x(i+1)- ...
                                  bubbleList{n}.x(i-1))* ...
                                 (fluidProp.dispRho-fluidProp.contRho); 
                    faceY = bubbleList{n}.distributeLagrangianToEulerian( ...
                            domain, faceY, bubbleList{n}.x(i), ...
                            bubbleList{n}.y(i), forceY, 2);
                end
            end
            
            % construct the density field using SOR
            for iter=1:param.maxIter
                oldRho = obj.rho;
                for i=2:domain.nx+1
                    for j=2:domain.ny+1
                        obj.rho(i,j) = (1.0-param.beta)*obj.rho(i,j)+ ...
                            param.beta*0.25* ...
                           (obj.rho(i+1,j  )+obj.rho(i-1,j  )+ ...
                            obj.rho(i  ,j+1)+obj.rho(i  ,j-1)+...
                            domain.dx*faceX(i-1,j  )- ...
                            domain.dx*faceX(i  ,j  )+ ...
                            domain.dy*faceY(i  ,j-1)- ...
                            domain.dy*faceY(i  ,j  ));
                    end
                end
                if max(max(abs(oldRho-obj.rho))) < param.maxErr
                    break
                end
            end              
        end
        
    end
end
