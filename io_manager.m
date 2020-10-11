% Input output manager.
classdef io_manager    
    methods(Static)
        % Clean output directory.
        function clean_dir
            delete output/bub_*;
        end
        
        % Create output directory.
        function create_dir
            mkdir output;
        end
    end
end
