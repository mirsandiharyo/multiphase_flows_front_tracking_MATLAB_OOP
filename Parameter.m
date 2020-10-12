% Parameter class.
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
        %% Initialize the simulation parameters.
        function obj = Parameter(nstep, dt, maxIter, maxErr, beta, outputFreq)
            obj.nstep = nstep;
            obj.dt = dt;
            obj.maxIter = maxIter;
            obj.maxErr = maxErr;
            obj.beta = beta;
            obj.outputFreq = outputFreq;
            obj.time = 0.0;
        end
        
        %% Increment the time.
        function incrementTime(obj)
            obj.time = obj.time+obj.dt;
        end
    end
end
