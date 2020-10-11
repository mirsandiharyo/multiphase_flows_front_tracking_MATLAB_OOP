% Parameter class.
classdef Parameter
    properties(SetAccess = private)
        nstep
        dt
        maxIter
        maxErr
        beta
        outFreq
    end
    
    methods
        % Initialize the simulation parameters.
        function obj = Parameter(nstep, dt, maxIter, maxErr, beta, outFreq)
            obj.nstep = nstep;
            obj.dt = dt;
            obj.maxIter = maxIter;
            obj.maxErr = maxErr;
            obj.beta = beta;
            obj.outFreq = outFreq;
        end
    end
end
