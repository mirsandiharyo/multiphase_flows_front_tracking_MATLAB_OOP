% Fluid properties class.
classdef FluidProp
    properties(SetAccess = private)
        rho
        mu
        sigma
    end
    
    methods
        % Initialize the simulation parameters.
        function obj = FluidProp(rho, mu, sigma)   
            obj.rho = rho;
            obj.mu = mu;
            obj.sigma = sigma;
        end
    end
end
