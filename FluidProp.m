% Fluid properties class.
classdef FluidProp
    properties(SetAccess = private)
        contRho
        contMu
        dispRho
        dispMu
        sigma
    end
    
    methods
        % Initialize the fluid properties of the continuous and 
        % dispersed phases.
        function obj = FluidProp(contRho, contMu, dispRho, dispMu, sigma)
            % Continuous phase
            obj.contRho = contRho;
            obj.contMu = contMu;
            % Dispersed phase            
            obj.dispRho = dispRho;
            obj.dispMu = dispMu;
            obj.sigma = sigma;
        end
    end
end
