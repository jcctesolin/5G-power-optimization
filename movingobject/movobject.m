% Parameters
numPeople = 50; % Number of people
numSteps = 100; % Number of steps in the simulation
areaSize = 100; % Size of the area (100x100)

% Initialize positions
positions = rand(numPeople, 2) * areaSize;

% Plot initial positions
figure;
scatter(positions(:,1), positions(:,2), 'filled');
xlim([0 areaSize]);
ylim([0 areaSize]);
title('People Movement Simulation');
hold on;

% Simulate movement
for step = 1:numSteps
    % Random movement: each person moves a small random step
    stepSize = 1; % Maximum step size
    movements = (rand(numPeople, 2) - 0.5) * stepSize * 2;
    positions = positions + movements;
    
    % Ensure positions stay within the area
    positions = max(min(positions, areaSize), 0);
    
    % Update plot
    scatter(positions(:,1), positions(:,2), 'filled');
    pause(0.1); % Pause to visualize movement
end

hold off;
