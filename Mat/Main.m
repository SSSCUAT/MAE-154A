%% 154A Project 
%known values 
% Conversion Facto
feet_meter = 1/3.28084;
% Aerodynamic Variable Definitions
CL = 1.5;    % 3D wing lift coefficient
AR = 7.4;    % Wing aspect ratio (Example value for Cessna 182)
CD_p = 0.027; %parasite Drag 
e = 0.80;            % Oswald's efficiency factor

% Induced Drag Calculation
% Formula from image: CDi = CL^2 / (pi * AR)
CDi = (CL^2) / (pi * AR);

%% Inputs - Change these values
v = 50*feet_meter;          % Velocity in meters per second (m/s)
L = 0.767*feet_meter;        % Chord length in meters (Suggested for Cessna 182)

%% Sea Level Constants (ISA Standard)
rho = 1.225;     % Air density at sea level (kg/m^3)
mu = 1.789e-5;   % Dynamic viscosity of air (kg/(m*s))

%% Calculation
% Re = (Density * Velocity * Length) / Viscosity
Re = (rho * v * L) / mu;

%% Display Output
fprintf('At Velocity: %.2f m/s\n', v);
fprintf('With Chord:  %.2f m\n', L);
fprintf('---------------------------\n');
fprintf('Reynolds Number: %.0f\n', Re);

% Display Result
fprintf('The Induced Drag Coefficient (CDi) is: %.4f\n', CDi);