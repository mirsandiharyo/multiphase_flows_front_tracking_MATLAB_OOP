% FluidProp class contains properties and methods related to the 
% fluid properties.
classdef FluidProp
    properties(SetAccess = private)
        contRho
        contMu
        dispRho
        dispMu
        sigma
    end
    
    methods
        %%
        function obj = FluidProp(contRho, contMu, dispRho, dispMu, sigma)
        % Initialize the fluid properties of the continuous and the
        % dispersed phases.    
            % Continuous phase.
            obj.contRho = contRho;
            obj.contMu = contMu;
            % Dispersed phase.           
            obj.dispRho = dispRho;
            obj.dispMu = dispMu;
            obj.sigma = sigma;
        end
    end
end
