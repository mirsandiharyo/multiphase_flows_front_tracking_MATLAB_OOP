% IOManager class contains static methods for input and output.
classdef IOManager    
    methods(Static)
        %% 
        function cleanDir
        % Clean output directory.    
            delete output/bub_*;
        end
        
        %%
        function createDir
        % Create output directory.    
            mkdir output;
        end
        
        %% 
        function [domain, param, fluidProp, bubbleList] = readInputFile
        % Read simulation parameters from the input file.    
            disp('choose the input file (.txt)');
            [inputName, filePath] = uigetfile('.txt');
            originPath = pwd;
            cd(filePath);
            fid = fopen(inputName);
            cd(originPath);
            
            % Solver parameters.
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
            
            % Numerical parameters.
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
            
            % Physical properties.
            % Dispersed phase.
            readLine = fgetl(fid);
            readLine = fgetl(fid);    
            readLine = regexp(fgetl(fid), '=', 'split');
            dispRho = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            dispMu = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            sigma = str2double(readLine{2});
            % Continuous phase.
            readLine = fgetl(fid);
            readLine = regexp(fgetl(fid), '=', 'split');
            contRho = str2double(readLine{2});
            readLine = regexp(fgetl(fid), '=', 'split');
            contMu = str2double(readLine{2});
            fluidProp = FluidProp(contRho, contMu, dispRho, dispMu, sigma);
            readLine = fgetl(fid);
            
            % Bubble size and location.
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
        
        %%
        function visualizeResults(domain, face, center, fluid, bubbleList, ...
            fluidProp, time, nstep)
        % Visualize the phase fraction field, velocity vector, velocity 
        % contour and marker points.
            % Calculate velocity at cell center.
            uCenter(1:domain.nx+1,1:domain.ny+1) = 0.5* ...
                (face.u(1:domain.nx+1,2:domain.ny+2)+ ...
                 face.u(1:domain.nx+1,1:domain.ny+1));
            vCenter(1:domain.nx+1,1:domain.ny+1) = 0.5* ...
                (face.v(2:domain.nx+2,1:domain.ny+1)+ ...
                 face.v(1:domain.nx+1,1:domain.ny+1));
            velMag = sqrt(uCenter.^2 + vCenter.^2);
            
            % Calculate phase fraction.    
            alpha = bsxfun(@minus,fluid.rho,fluidProp.contRho);
            alpha = bsxfun(@times,alpha,1/ ...
                          (fluidProp.dispRho-fluidProp.contRho));
            
            % Create the grid.
            gridX = linspace(0, domain.nx, domain.nx+1)*domain.dx;
            gridY = linspace(0, domain.ny, domain.ny+1)*domain.dy;
            hold off, 
            
            % Plot contour of velocity magnitude. 
            contour(gridX,gridY,flipud(rot90(velMag)),'linewidth',1.5), ...
            axis equal, axis([0 domain.lx 0 domain.ly]);
            hold on;
            
            % Plot the phase field.
            imagesc(center.x,center.y,flipud(rot90(alpha)),'AlphaData',0.9),
            colormap('jet'),colorbar,caxis([0 1]);
            
            % Set colorbar title.
            ph = colorbar;
            colorTitleHandle = get(ph,'Title');
            caption = 'Phase fraction';
            set(colorTitleHandle ,'String',caption,'FontSize', 10);
            
            % Plot velocity vector.
            quiver(gridX,gridY,flipud(rot90(uCenter)), ...
                flipud(rot90(vCenter)),'w');
            
            % Plot the marker points.
            for n=1:length(bubbleList)
                h = plot(bubbleList{n}.x(1:bubbleList{n}.point), ...
                         bubbleList{n}.y(1:bubbleList{n}.point), ...
                         'k','linewidth',2);
            end
            
            % Set title.
            caption = sprintf('Time = %f s', time);
            title(caption, 'FontSize', 10);     
            pause(0.001)
            
            % Save the plot.  
            caption = sprintf('output/bub_%03d.png',nstep);
            saveas(h, caption);
        end

    end
end
