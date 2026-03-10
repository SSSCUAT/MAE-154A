clear; clc; close all;

% ATMOSPHERE PROPERTIES ;
inputs.rho = 0.0023769;    % Air density [slug/ft^3] at sea level
inputs.mu  = 3.737e-7;     % Dynamic viscosity of air [slug/(ft*s)]
inputs.a   = 1116.45;      % Speed of sound [ft/s] at sea level
% BASE AIRCRAFT GEOMETRY & AERODYNAMIC VALUES
% Wing geometry
inputs.Arw     = 7.4;        
inputs.bw      = 5.617;      
inputs.lambdaw = 0.40;       

% Horizontal tail geometry
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
inputs.EF    = 0.7;              
inputs.W     = 23;               

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
inputs.V      = linspace(30, 300, 1000); 

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
inputs.Wfuel       = 0.5;   %Pounds  
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
inputs.x_fuel      = 1.25
% Endurance 
inputs.c_p_hp_hr = 2.2;  % Specific fuel consumption [lb/(hp*hr)]

%% =====================================
% DISPLAY / OUTPUT OPTIONS
% =====================================
inputs.makePlots            = 0;   
inputs.makeTable            = 0;   
inputs.makePrint_stability  = 0;   
inputs.makePrint_tail_Volume = 0;


%% =====================================
% MONTE CARLO DESIGN STUDY
% =====================================

% Number of random aircraft designs to test
nCases = 10000;

%% VARIABLE RANGES (EDIT THESE WHENEVER YOU WANT TO CHANGE DESIGN SPACE)

% (Extracted from 5-design set)
L_fuse_min = 3.5187;   L_fuse_max = 3.6212;
Arw_min    = 6.3898;   Arw_max    = 6.5140;
Art_min    = 3.1548;   Art_max    = 3.4484;
bw_min     = 4.5610;   bw_max     = 4.6906;
x_wle_min  = 1.0540;   x_wle_max  = 1.0841;


%% STORAGE ARRAYS
L_fuse_vals = zeros(nCases,1);
Arw_vals    = zeros(nCases,1);
Art_vals    = zeros(nCases,1);
bw_vals     = zeros(nCases,1);
x_wle_vals  = zeros(nCases,1);

Drag        = zeros(nCases,1);
StallSpeed  = zeros(nCases,1);
MaxSpeed    = zeros(nCases,1);
ROC_Stall   = zeros(nCases,1);
W_total     = zeros(nCases,1);
SM_dry      = zeros(nCases,1);
SM_0        = zeros(nCases,1);
Endurance_min = zeros(nCases,1);

% Performance constraints
min_StallSpeed = 55;   % must be BELOW this
min_MaxSpeed   = 120;  % must be ABOVE this
min_ROC        = 40;   % must be ABOVE this

%% MONTE CARLO LOOP
% IMPORTANT: Make sure your baseline 'inputs' struct is loaded/defined 
% BEFORE this loop if computeAircraft needs other constant values!

for i = 1:nCases
    
    % HARDCODED IDEAL AIRCRAFT (Comment out the rand equations)
    inputs.L_fuse = L_fuse_min + rand*(L_fuse_max - L_fuse_min);
    inputs.Arw    =  Arw_min    + rand*(Arw_max - Arw_min);
    inputs.Art    = Art_min    + rand*(Art_max - Art_min);
    inputs.bw     = bw_min     + rand*(bw_max - bw_min);
    inputs.x_wle  =   x_wle_min  + rand*(x_wle_max - x_wle_min);

    % Store the generated design variables
    L_fuse_vals(i) = inputs.L_fuse;
    Arw_vals(i)    = inputs.Arw;
    Art_vals(i)    = inputs.Art;
    bw_vals(i)     = inputs.bw;
    x_wle_vals(i)  = inputs.x_wle;

    % Run aircraft model
    outputs = computeAircraft(inputs);

    % Store outputs
    StallSpeed(i) = outputs.V_stall;
    MaxSpeed(i)   = outputs.V_max;
    ROC_Stall(i)  = outputs.ROC_stall;
    SM_0(i)       = outputs.SM_0;
    SM_dry(i)     = outputs.SM_dry;
    W_total(i)    = outputs.W_total;
    Endurance_min(i) = outputs.Endurance_min;
end

%% ==============================
% SAVE ALL CASES INTO ONE CSV FILE
% ==============================

stall_ok = StallSpeed <= min_StallSpeed;
max_ok   = MaxSpeed   >= min_MaxSpeed;
roc_ok   = ROC_Stall  >= min_ROC;

SM_0_ok   = SM_0 >= 0.04 & SM_0 <= 0.18;
SM_dry_ok = SM_dry >= 0.04 & SM_dry <= 0.18;

E_ok = Endurance_min >= 45;

PassFlag = stall_ok & max_ok & roc_ok & SM_dry_ok & SM_0_ok & E_ok;

%% TEXT CLEANER
Result = strings(length(PassFlag),1);
Result(PassFlag)  = "PASS";
Result(~PassFlag) = "FAIL";

folderName = fullfile(pwd,'Results');  

% Safely create the Results folder if it doesn't exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

%% CREATE TABLE (DESIGN VARIABLES FIRST)
ResultsTable = table( ...
    L_fuse_vals, x_wle_vals, Arw_vals, Art_vals, bw_vals, SM_0, SM_dry, ...
    MaxSpeed, ROC_Stall, StallSpeed, Endurance_min, W_total, Result);

% Assign the column headers in the exact same order
ResultsTable.Properties.VariableNames = ...
    {'L_fuse','x_wle', 'Arw', 'Art', 'bw', 'SM_0', 'SM_dry', ...
    'Max_Speed', 'ROC_Stall', 'Stall_Speed', 'Endurance_min', 'W_total', 'Result'};

%% CREATE TIMESTAMP
timestamp = string(datetime("now", 'Format', 'yyyy-MM-dd_HH-mm-ss'));

%% FILE NAME
fileName = fullfile(folderName, ...
    sprintf('MonteCarlo_%s.csv', timestamp));

%% SAVE CSV
writetable(ResultsTable,fileName);

%% DISPLAY CONFIRMATION
disp(['All cases saved to: ' fileName])