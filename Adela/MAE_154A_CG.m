clear, clc;

run('MAE_154A_weight.m');

%% Plane specs

rho = 0.0023769;  % air density (slug/ft^3)
V_stall = 50;        % stall velocity (ft/s)
CL_MAX = 1.5;      % Cl max 
%W = 19;        % total weight of aircraft + avioinics/payload/etc. (lbs)
W = W_pick;

%% Plane Parameter

 L_fuse = 5; % length of fuselage (ft)

%% Wing Parameters

Arw = 7.4;         % aspect ratio
Sw = 4.26; % wing area
bw = 5.6; % wing span (tip to tip)
c = Sw/bw;          % chord length 
%c = 0.8; 

%% Horizontal tail parameters
lambdaw =0.4;
lambdath = 0.5;
Sth = 0.2 * Sw;
Art = 4;          % horizontal tail aspect ratio (GUESSED!!)
bth = sqrt(Art*Sth); % Horizontal tail wing span
V_ht = 0.5;       % horizontal tail volume coefficient (SEARCHED UP!!)

Crw = (2*Sw)/(bw*(1 + lambdaw));
Ctw = lambdaw * Crw;
Cw  = (2/3) * Crw * ((1 + lambdaw + lambdaw^2)/(1 + lambdaw)); % wing MAC

Crth = (2*Sth)/(bth*(1 + lambdath));

x_wle = 1;

hacw = 0.25;   % wing AC at 0.25 MAC
hact = 0.25;   % tail AC at 0.25 chord (root chord ref)

x_wac = x_wle + hacw*Cw;                 % wing AC from nose [ft]
x_tac = L_fuse - (1 - hact)*Crth;        % tail AC from nose [ft]
L_ht   = x_tac - x_wac;                   % wing AC -> tail AC [ft]

% L_ht = 4*c;         % distance from AC of wing to AC of horizontal tail
c_ht = 0.38;
% c_ht = Sth / bth;     % Horizontal tail chord length
% Sth = (V_ht*Sw*c) / L_ht; % Horizonatil tail area

%% Vertical tail parameters ( Sv, bv, cv calculated further down)

AR_vt = 1.8;      % vertical tail aspect ratio (CONSTANT!!)
V_vt = 0.04;      % volume coefficient vertical tail (CONSTANT!!)
% Temporary!!!
Lv = 4*c;         % distance from AC of wing to AC of verticla tail

%% Structural Weights (use from weight code)

Wstruct = 4.3720; % total weight of empty aircrafts (lbs)
Wf = 0.6465;      % weight of fuselage (lbs)
Ww = 1.2398;      % Weight of wing
Wht = 0.6363;  % weight of horizontal tail
Wvt = 0.3493;  % weight of vertical tail
Wlg = 1.5;     % weight of landing gear
Wprop = .24;     % weight of propellar

%% locations of LE, AC, and CG of wing, htail, vtail, landing gear, & prop

x_wle = 1; % wing leading edge wrt nose (ft)
x_wac = x_wle + 0.25*c; % location of wing aerodynamic center wrt nose
x_wcg = x_wle + 0.35*c; % location of wing center of gravity wrt nose

x_hac = x_wac + L_ht;        % location of horizontal tail AC wrt nose
x_hle = x_hac - 0.25*c_ht; % location of horizontal tail leading edge wrt nose
x_hcg = x_hle + 0.35*c_ht; % location of horizontal tail cg wrt nose

x_vac = x_wac + Lv;         % location of vertical tail AC wrt nose

x_fcg = 0.4*L_fuse; % CG of fuselage (estimated)
 x_lgcg = 2; % CG of landing gear (estimated)
 x_propcg = -0.2; % CG of propellar

%% Vertical Tail (contd)
W_novt = [Ww, Wht, Wf, Wlg, Wprop];  % weight without vertical tail
x_novt = [x_wcg, x_hcg, x_fcg, x_lgcg, x_propcg];  

x_cg_noVT = sum(W_novt.*x_novt) / sum(W_novt); % cg without vertical tail

% Vertical Tail Parameters
l_vt = x_vac - x_cg_noVT % distance from cg wihtout vt to ac of vt 
S_vt = V_vt*((Sw*bw)/l_vt) % VERTICAL TAIL AREA
b_vt = sqrt(AR_vt*S_vt)    % VERTICAL TAIL span
c_vt = S_vt / b_vt           % VERTICAL TAIL chord (NACA 0012)

x_vle = x_vac - 0.25*c_vt;  % location of vertical tail leading edge wrt nose
x_vcg = x_vle + 0.35*c_vt; % location of vertical tail cg

%% Total Center of gravity of structure (empty plane) !!!
W = [Ww, Wht, Wvt, Wf, Wlg, Wprop];
x = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_propcg];

x_cg_str = sum(W.*x) / sum(W); % Total structural center of gravity 

%% Weight of Avionics (lbs)
Wcam = 0.95;  % weight of camera 
Wcomp = 0.10; % weight of computer
Wgps = 0.07;  % weight of GPS
Wbat = 0.44;  % weight of battery
Wserv = 0.06; % weight of servos
Weng = 1.76;     % weight of engine
Wfs = 0.25;   % weight of fuel system

%% Weight of payload/fuel
Wbal = 0.66;   % weight of 1 waterballoon
N_bal = 3;
Wpay = N_bal*Wbal; % weight of payload
%Wpay = 0;
Wfuel = 1.5;   % weight of initial fuel
Wfuel_half = 0.5*Wfuel;
Wfuel_quart = 0.25*Wfuel;

%% Locations of avionics wrt to nose
x_cam = 0.78;    % location of camera
x_comp = 2.5; % location of computer
x_gps = 2.5;  % location of GPS
x_bat = 1;  % location of battery
x_serv = 1.25; % location of servos
x_eng = 0.25; % location of engine
x_fs = 1.25;   % location of fuel system 

%% Locations of payload/fuel wrt nose
x_1bal = 1.25;
X_2bal = 1.25;
x_3bal = 1.25;
x_pay = 1.25;  % location of payload
x_fuel = 1.25; % loaction of fuel 

%% Total CG without payload or fuel!!
Wdry = [W, Wcam, Wcomp, Wgps, Wbat, Wserv, Weng, Wfs]; %weight of structure + avionics
x_dry  = [x, x_cam, x_comp, x_gps, x_bat, x_serv, x_eng, x_fs];

x_cg_dry = sum(Wdry.*x_dry)/sum(Wdry)

%% Initial CG with payload and fuel!!
W0 = [Wdry, Wpay, Wfuel]; % initial weight w/ payload and fuel
x_0  = [x_dry, x_pay, x_fuel];

x_cg_0 = sum(W0.*x_0)/sum(W0) % initial cg w/ payload and fuel

%% CG with 1/2 fuel and all payload
W_fhalf_pay = [Wdry, Wpay, Wfuel_half]; % weight w/ payload and 1/2 fuel
x_fhalf_pay  = [x_dry, x_pay, x_fuel];

x_cg_fhalf_pay = sum(W_fhalf_pay.*x_fhalf_pay)/sum(W_fhalf_pay) % cg w/ payload and 1/2 fuel

%% CG with 1/2 fuel and no payload
W_fhalf_nopay = [Wdry, Wfuel_half]; % weight w/ no payload and 1/2 fuel
x_fhalf_nopay  = [x_dry, x_fuel];

x_cg_fhalf_nopay = sum(W_fhalf_nopay.*x_fhalf_nopay)/sum(W_fhalf_nopay) % cg w/ no payload and 1/2 fuel

%% CG with 1/4 fuel and no payload
W_fquart_nopay = [Wdry, Wfuel_quart]; % weight w/ no payload and 1/4 fuel
x_fquart_nopay  = [x_dry, x_fuel];

x_cg_fquart_nopay = sum(W_fquart_nopay.*x_fquart_nopay)/sum(W_fquart_nopay) % cg w/ no payload and 1/4 fuel

%% CG with fuel and no payload
W_f_nopay = [Wdry, Wfuel]; % weight w/ no payload and all fuel
x_f_nopay  = [x_dry, x_fuel];

x_cg_f_nopay = sum(W_f_nopay.*x_f_nopay)/sum(W_f_nopay) % cg w/ no payload and all fuel
