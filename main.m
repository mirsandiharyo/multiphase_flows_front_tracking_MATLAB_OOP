%% 
% A two-dimensional gas-liquid multiphase flows using a front-tracking type
% method. A set of Navier-Stokes equation is solved on a eulerian grid 
% using a second order projection method. The fluid properties are advected 
% by lagrangian marker points. The time marching is second order by using 
% predictor-corrector method. The code can be used to simulate a bubble 
% rising in a rectangular box.
% Created by: Haryo Mirsandi
clear;
close all;
%% Initialization
% Clean output folder
IOManager.createDir;
IOManager.cleanDir;

% read input file
[domain, param, fluidProp, bubbleList] = IOManager.readInputFile();

% initialize variables (grid, velocity, pressure, and force)
face = Face(domain);
center = Center(domain);

% initialize the physical properties inside the domain
fluid = Fluid(domain, fluidProp);
fluid.initializeDomain(domain, center, bubbleList, fluidProp);

% set the initial front (gas-liquid interface)
for n=1:length(bubbleList)
    bubbleList{n}.initializeFront();
end
                   
% visualize the initial condition
IOManager.visualizeResults(domain, face, center, fluid, bubbleList, ...
            fluidProp, param.time, 0)

%% start time loop
for nstep=1:param.nstep
    
    % store second order variables
    face.storeOldVariables();
    fluid.storeOldVariables();
    for n=1:length(bubbleList)
        bubbleList{n}.storeOldVariables();
    end

        % calculate the surface tension force at the front (lagrangian grid)
        % and distribute it to eulerian grid

        % update the tangential velocity at boundaries

        % calculate the (temporary) velocity

        % solve pressure
        
        % correct the velocity by adding the pressure gradient
        
        % update the front location 

        % update physical properties
  
    % end
    
    % store second order variables
    face.store2ndOrderVariables();
    fluid.store2ndOrderVariables();
    for n=1:length(bubbleList)
        bubbleList{n}.store2ndOrderVariables();
    end 
    
    % restructure the front
    for n=1:length(bubbleList)
        bubbleList{n}.restructure_front(domain);
    end     

    % visualize the results
    param.incrementTime();
    if mod(nstep, param.outputFreq) == 0
        IOManager.visualizeResults(domain, face, center, fluid, bubbleList, ...
            fluidProp, param.time, nstep)
    end
end
%% end time loop
disp('program finished');

