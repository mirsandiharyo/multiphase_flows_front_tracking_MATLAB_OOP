% Parameter class contains properties and methods related to the 
% computational parameters.
classdef Parameter < handle
    properties(SetAccess = private)
        nstep
        dt
        maxIter
        maxErr
        beta
        outputFreq
        time
    end
    
    methods
        %% 
        function obj = Parameter(nstep, dt, maxIter, maxErr, beta, ...
                outputFreq)
        % Initialize the simulation parameters.    
            obj.nstep = nstep;
            obj.dt = dt;
            obj.maxIter = maxIter;
            obj.maxErr = maxErr;
            obj.beta = beta;
            obj.outputFreq = outputFreq;
            obj.time = 0.0;
        end
        
        %% 
        function incrementTime(obj)
        % Increment the time.    
            obj.time = obj.time+obj.dt;
        end
    end
end
