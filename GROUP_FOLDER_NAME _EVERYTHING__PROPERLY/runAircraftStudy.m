clear; clc; close all;

%% =====================================
% ATMOSPHERE PROPERTIES
% =====================================
inputs.rho = 0.0023769;    % Air density [slug/ft^3] at sea level
inputs.mu  = 3.737e-7;     % Dynamic viscosity of air [slug/(ft*s)]
inputs.a   = 1116.45;      % Speed of sound [ft/s] at sea level

%% =====================================
% BASE AIRCRAFT GEOMETRY & AERODYNAMIC VALUES
% =====================================
% Wing geometry
inputs.Arw     = 7.4;        
inputs.bw      = 5.617;      
inputs.Sw      = 4.264931025; 
inputs.lambdaw = 0.40;       

% Horizontal tail geometry
inputs.Sth      = 0.5331163781; 
inputs.bth      = 1.460296378;  
inputs.lambdath = 0.50;         
inputs.Art      = 4.0;          

% Aerodynamic parameters
inputs.e        = 0.85;         
inputs.cla      = 5.73;         
inputs.downwash = 0.25;         
inputs.it       = -1*pi/180;    
inputs.Cmacw    = -0.05;        
inputs.tau      = 0.5;          
inputs.CL_MAX   = 1.5;          

% Fuselage geometry
inputs.L_fuse = 5.0;             
inputs.W_fuse = 8/12;            
inputs.H_fuse = 8/12; 
inputs.D_fuse = 1/24;

% Propulsion / engine
inputs.Power = 2.8;              
inputs.EF    = 0.6;              
inputs.W     = 20;               

% NOTE: Base x_wle (0.75) and hardcoded CGs were removed here because 
% Adela's section updates x_wle to 1.0, and the function calculates the CGs dynamically!

%% =====================================
% JONATHAN'S VARIABLES
% =====================================
% NOTE: Jonathan's hardcoded CGs and AR_vt (1.8) were removed so they don't 
% conflict with Adela's dynamic calculations and updated AR_vt (1.5).

% new vertical tail components 
inputs.c_vt   = 0.5;            
inputs.tc_v   = 0.12;            
inputs.xcm_v  = 0.30;            
inputs.Lam_v  = 0;               
inputs.V      = linspace(30, 100, 500); 

%% =====================================
% ADELA'S INPUTS
% =====================================

% Horizonatl Tail 
inputs.c_ht        = 0.38;
inputs.hac_ht      = 0.25;   % x_ac/c
inputs.V_ht        = 0.5;    % ht volume coeff

% Vertical Tail
inputs.hacv        = 0.25;
inputs.AR_vt       = 1.8;
inputs.V_vt        = 0.04;   % vt volume coeff
inputs.Sv          = 0.4;
inputs.bv          = 0.7;

% Wing 
inputs.c           = 0.8;
inputs.hac                = 0.25;

% for Weight code
inputs.Wguess      = 25;      
inputs.V_stall     = 50; 
inputs.N_load      = 6.6;
inputs.tc          = 0.12;
inputs.V_max_kts   = 80;           

% Specific Component Weights [lbs]
inputs.Wlg         = 1.5;
inputs.Wprop       = 0.24;
inputs.Weng        = 1.76;
inputs.Neng        = 1;       
inputs.Wfuel       = 1.5;     
inputs.Wfs         = 0.25;
inputs.Wau         = 1.62;
inputs.Wbal        = 0.66;
inputs.N_bal       = 3;
inputs.Wcam        = 0.95;    
inputs.Wcomp       = 0.10;    
inputs.Wgps        = 0.07;    
inputs.Wbat        = 0.44;    
inputs.Wserv       = 0.06;           

% Component CG Locations (wrt nose) [ft]
inputs.x_wle       = 1.0;     % This overrides the old base value
inputs.x_fcg       = 0.4 * inputs.L_fuse; 
inputs.x_lgcg      = 2.0;     
inputs.x_propcg    = -0.2;    
inputs.x_cam       = 0.78;    
inputs.x_comp      = 2.5;     
inputs.x_gps       = 2.5;     
inputs.x_bat       = 1.0;     
inputs.x_serv      = 1.25;    
inputs.x_eng       = 0.25;    
inputs.x_fs        = 1.25;    
inputs.x_pay       = 1.25;    
inputs.x_fuel      = 1.25;

%% =====================================
% DISPLAY / OUTPUT OPTIONS
% =====================================
inputs.makePlots            = 0;   
inputs.makeTable            = 0;   
inputs.makePrint_stability  = 0;   
inputs.makePrint_tail_Volume = 0;

%% VARIABLE SWEEP

sweep_var_name = 'e';               % Name of the field in 'inputs' to vary (e.g., 'W', 'Power', 'EF')
sweep_values   = [0.5:0.01:1];     % Values to sweep through for that variable
%% STORAGE ARRAYS (ADD THIS ABOVE YOUR SWEEP LOOP)

nCases = length(sweep_values);

SweepValue  = zeros(nCases,1);
Drag      = zeros(nCases,1);
StallSpeed  = zeros(nCases,1);
MaxSpeed    = zeros(nCases,1);
ROC_Stall   = zeros(nCases,1);
%
min_StallSpeed = 55;   % must be BELOW this
min_MaxSpeed   = 120;   % must be ABOVE this
min_ROC        = 67;    % must be ABOVE this
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
    % Store results for this case
SweepValue(i) = sweep_values(i);
StallSpeed(i) = outputs.V_stall;
MaxSpeed(i)   = outputs.V_max;
ROC_Stall(i)  = outputs.ROC_stall;
end
%% ==============================
% SAVE ALL SWEEP CASES INTO ONE CSV FILE
% ==============================
% Logical checks (returns true/false arrays)
stall_ok = StallSpeed <= min_StallSpeed;
max_ok   = MaxSpeed   >= min_MaxSpeed;
roc_ok   = ROC_Stall  >= min_ROC;
PassFlag = stall_ok & max_ok & roc_ok; 
%%text cleaner
Result = strings(length(PassFlag),1);

Result(PassFlag) = "PASS";
Result(~PassFlag) = "FAIL";

% fullfile(pwd,'Results')
% - pwd gives the current working directory (where your script is running)
% - fullfile safely builds a path using correct slashes for your OS
% Result: ./Results

folderName = fullfile(pwd,'Results');  


% Create a table object
% Each input column must be the same length
% Since SweepValue, StallSpeed, etc. were filled inside the loop,
% each row corresponds to ONE sweep case
ResultsTable = table(SweepValue, StallSpeed, MaxSpeed, ROC_Stall,Result);


% Rename the table column headers
% The first column name is dynamic (whatever you set sweep_var_name to)
% Example: if sweep_var_name = 'Arw', first column becomes 'Arw'
ResultsTable.Properties.VariableNames = ...
    {sweep_var_name,'Stall_Speed','Max_Speed','ROC_Stall','Result'};


% Create a timestamp string so each run creates a NEW file
% now         → current date/time
% datestr     → converts it to readable string
% Format: yyyy-mm-dd_HH-MM-SS
% Example: 2026-02-23_14-42-11
timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');


% Build the full filename
% sprintf builds a formatted string:
%   Simulation_<SweepVariable>_<Timestamp>.csv
%
% Example output:
%   Simulation_Arw_2026-02-23_14-42-11.csv
fileName = fullfile(folderName, ...
    sprintf('Simulation_%s_%s.csv', sweep_var_name, timestamp));


% Write the table to a CSV file
% writetable automatically:
%   - Writes column headers
%   - Writes each row
%   - Handles formatting
writetable(ResultsTable,fileName);


% Display confirmation in command window
% Helps you verify where it saved
disp(['All cases saved to: ' fileName])