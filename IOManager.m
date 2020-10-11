% Input output manager.
classdef IOManager    
    methods(Static)
        % Clean output directory.
        function cleanDir
            delete output/bub_*;
        end
        
        % Create output directory.
        function createDir
            mkdir output;
        end
    end
end
