clear; clc; close all;

%% =====================================
% ATMOSPHERE PROPERTIES
% =====================================
% These are the atmospheric parameters used to compute aerodynamic forces
inputs.rho = 0.0023769;    % Air density [slug/ft^3] at sea level (standard day)
inputs.mu  = 3.737e-7;     % Dynamic viscosity of air [slug/(ft*s)]
inputs.a   = 1116.45;      % Speed of sound [ft/s] at sea level (used for Mach number)

%% =====================================
% BASE AIRCRAFT GEOMETRY & AERODYNAMIC VALUES
% =====================================
% Wing geometry
inputs.Arw     = 7.4;        % Aspect ratio of wing (span^2 / wing area)
inputs.bw      = 5.617;      % Wing span [ft]
inputs.Sw      = 4.264931025; % Wing planform area [ft^2]
inputs.lambdaw = 0.40;       % Wing taper ratio (tip chord / root chord)

% Horizontal tail geometry
inputs.Sth      = 0.5331163781; % Tail planform area [ft^2]
inputs.bth      = 1.460296378;  % Tail span [ft]
inputs.lambdath = 0.50;         % Tail taper ratio
inputs.Art      = 4.0;          % Tail aspect ratio

% Aerodynamic parameters
inputs.e        = 0.85;         % Oswald efficiency factor (for induced drag)
inputs.cla      = 5.73;         % 2D lift curve slope [per radian]
inputs.downwash = 0.25;         % Downwash effect from wing at tail
inputs.it       = -1*pi/180;    % Tail incidence angle [radians] negative for nose-down
inputs.Cmacw    = -0.05;        % Wing pitching moment coefficient about aerodynamic center
inputs.tau      = 0.5;          % Elevator effectiveness factor
inputs.CL_MAX   = 1.5;          % Maximum lift coefficient (for stall calculation)

% Longitudinal reference locations
inputs.x_wle    = 0.75;         % Wing leading edge position along fuselage [ft]
inputs.x_cg_str = 1.7147;       % CG location with structural load [ft]
inputs.x_cg_tot = 1.4103;       % Total CG location including payload/fuel [ft]
inputs.x_cg_0   = 1.5071;       % Initial CG location [ft]

% Fuselage geometry
inputs.L_fuse = 5.0;             % Fuselage length [ft]
inputs.W_fuse = 8/12;            % Fuselage width [ft]
inputs.H_fuse = 8/12;            % Fuselage height [ft]

% Propulsion / engine
inputs.Power = 2.8;              % Available power [hp]
inputs.EF    = 0.6;              % Propeller efficiency (0–1)

% Default aircraft weight (will be overwritten by sweep)
inputs.W = 20;                    % Aircraft weight [lb]

%% =====================================
% DISPLAY / OUTPUT OPTIONS
% =====================================
% Control what the computeAircraft function plots or prints
inputs.makePlots            = 0;   % 1 = create figures, 0 = skip plots
inputs.makeTable            = 0;   % 1 = display table of results
inputs.makePrint_stability  = 0;   % 1 = print static margin and stability info
inputs.makePrint_tail_Volume = 0;  % 1 = print tail volume info

%% VARIABLE SWEEP

sweep_var_name = 'W';               % Name of the field in 'inputs' to vary (e.g., 'W', 'Power', 'EF')
sweep_values   = [22 25 28 30];     % Values to sweep through for that variable

% Loop over each value in the sweep_values array
for i = 1:length(sweep_values)
    
    % Dynamically update the chosen input variable in the inputs struct
    inputs.(sweep_var_name) = sweep_values(i);  

    % Call the main computational function that computes aerodynamics,
    % stability, drag, lift, power required, ROC, stall speed, max speed, etc.
    outputs = computeAircraft(inputs);

    % Display results in the command window for each sweep value
    % This prints: sweep variable, max L/D velocity, stall speed, max speed,
    % and rate of climb at stall speed
    fprintf("%s = %.2f | V_LDmax = %.2f | V_stall = %.2f | V_max = %.2f | ROC_stall = %.2f ft/s\n", ...
        sweep_var_name, sweep_values(i), outputs.V_LDmax, outputs.V_stall, outputs.V_max, outputs.ROC_stall);

    % Optional: plot Rate of Climb vs Velocity for this sweep value
    % Uncomment to see ROC curves
    % figure;
    % plot(outputs.V, outputs.ROC, 'LineWidth', 2);
    % xlabel('Velocity [ft/s]'); ylabel('Rate of Climb [ft/s]');
    % title(sprintf('Rate of Climb vs Velocity (%s = %.2f)', sweep_var_name, sweep_values(i)));
    % grid on;

end