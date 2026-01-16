%% 154A Project 
%known values 

% Aerodynamic Variable Definitions
CL = 1.5;    % 3D wing lift coefficient
AR = 7.4;    % Wing aspect ratio (Example value for Cessna 182)
CD_p = 0.027; %parasite Drag 
e = 0.80;            % Oswald's efficiency factor

% Induced Drag Calculation
% Formula from image: CDi = CL^2 / (pi * AR)
CDi = (CL^2) / (pi * AR);

% Display Result
fprintf('The Induced Drag Coefficient (CDi) is: %.4f\n', CDi);