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
        
        % Read simulation parameters from the input file
        function [domain, param, fluidProp, bubbleList] = readInputFile
            disp('choose the input file (.txt)');
            [inputName, filePath] = uigetfile('.txt');
            originPath = pwd;
            cd(filePath);
            fid = fopen(inputName);
            cd(originPath);
            
            % Solver parameters
            readLine = fgetl(fid); %#ok<*NASGU>
            readLine = regexp(fgetl(fid), '=', 'split');
            nstep = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            dt = str2double(readLine{2}); 
            readLine = regexp(fgetl(fid), '=', 'split');
            maxIter = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            maxErr = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            beta = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            outputFreq = str2double(readLine{2});
            param = Parameter(nstep, dt, maxIter, maxErr, beta, outputFreq);
            readLine = fgetl(fid);
            readLine = fgetl(fid);
            
            % Numerical parameters
            readLine = regexp(fgetl(fid), '=', 'split');
            lx = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            ly = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            nx = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            ny = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            gravx = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            gravy = str2double(readLine{2});
            domain = Domain(lx, ly, nx, ny, gravx, gravy);
            readLine = fgetl(fid);
            
            % Physical properties
            % Dispersed phase
            readLine = fgetl(fid);
            readLine = fgetl(fid);    
            readLine = regexp(fgetl(fid), '=', 'split');
            dispRho = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            dispMu = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            sigma = str2double(readLine{2});
            % Continuous phase
            readLine = fgetl(fid);
            readLine = regexp(fgetl(fid), '=', 'split');
            contRho = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            contMu = str2double(readLine{2});
            fluidProp = FluidProp(contRho, contMu, dispRho, dispMu, sigma);
            readLine = fgetl(fid);
            
            % Bubble size and location
            readLine = fgetl(fid);
            readLine = regexp(fgetl(fid), '=', 'split');
            numBubble = str2double(readLine{2});
            bubbleList = cell(numBubble, 1);
            for i=1:numBubble
                readLine = regexp(fgetl(fid), '=', 'split');
                radius = str2double(readLine{2});
                readLine = regexp(fgetl(fid), '=', 'split');
                centerX = str2double(readLine{2});
                readLine = regexp(fgetl(fid), '=', 'split');
                centerY = str2double(readLine{2});
                readLine = regexp(fgetl(fid), '=', 'split');
                point = str2double(readLine{2});
                bubbleList{i} = Bubble(centerX, centerY, radius, point);
            end
            fclose(fid);
        end
        
    end
end
