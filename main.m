%% 
% A two-dimensional gas-liquid multiphase flows using a front-tracking type
% method. A set of Navier-Stokes equation is solved on a eulerian grid 
% using a second order projection method. The fluid properties are advected 
% by lagrangian marker points. The time marching is second order by using 
% predictor-corrector method. The code can be used to simulate a bubble 
% rising in a rectangular box.
% Created by: Haryo Mirsandi

%% Initialization
% Clean output folder.
IOManager.createDir;
IOManager.cleanDir;

% Read input file.
[domain, param, fluidProp, bubbleList] = IOManager.readInputFile();

% Initialize variables (grid, velocity, pressure, and force).
face = Face(domain);
center = Center(domain);

% Initialize the physical properties inside the domain.
fluid = Fluid(domain, fluidProp);
fluid.initializeDomain(domain, center, bubbleList, fluidProp);

% Set the initial front (gas-liquid interface).
for n=1:length(bubbleList)
    bubbleList{n}.initializeFront();
end
                   
% Visualize the initial condition.
IOManager.visualizeResults(domain, face, center, fluid, bubbleList, ...
            fluidProp, param.time, 0)

%% Start time loop
for nstep=1:param.nstep
    
    % Store second order variables.
    face.storeOldVariables();
    fluid.storeOldVariables();
    for n=1:length(bubbleList)
        bubbleList{n}.storeOldVariables();
    end

    % Second order loop.
    for substep=1:2
        
        % Calculate the surface tension force at the front (lagrangian grid)
        % and distribute it to eulerian grid.
        face.initializeForce(domain);   
        for n=1:length(bubbleList)
            bubbleList{n}.calculateSurfaceTension(domain, fluidProp, face);
        end 
            
        % Update the tangential velocity at boundaries.
        face.updateWallVelocity(domain);
        
        % Calculate the (temporary) velocity.
        face.calculateTemporaryVelocity(param, domain, fluidProp, fluid);
        
        % Solve the pressure field.
        center.solvePressure(domain, param, fluid, face);
        
        % Correct the velocity field by adding the pressure gradient.
        face.correctVelocity(domain, param, center, fluid);
        
        % Update the front location.
        for n=1:length(bubbleList)
            bubbleList{n}.updateFrontLocation(face, param, domain);
        end 
        
        % Update the physical properties.
        fluid.updateDensity(param, domain, bubbleList, fluidProp);
        fluid.updateViscosity(fluidProp);
    end
    
    % Store second order variables.
    face.store2ndOrderVariables();
    fluid.store2ndOrderVariables();
    for n=1:length(bubbleList)
        bubbleList{n}.store2ndOrderVariables();
    end 
    
    % Restructure the front.
    for n=1:length(bubbleList)
        bubbleList{n}.restructureFront(domain);
    end     

    % Visualize the results.
    param.incrementTime();
    if mod(nstep, param.outputFreq) == 0
        IOManager.visualizeResults(domain, face, center, fluid, bubbleList, ...
            fluidProp, param.time, nstep)
    end
end
%% end time loop
disp('program finished');

