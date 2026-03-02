clear, clc;

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
% Plane specs
r = 0.002377;  % air density (slug/ft^3)
V = 55;        % stall velocity (ft/s)
Cl = 1.5;      % Cl max 
W = 23;        % total weight of aircraft + avioinics/payload/etc. (lbs)

% Plane Parameter
Lf = 4; % length of fuselage (ft)

% Wing Parameters
AR = 7.4;         % aspect ratio
S = (W*2) / (Cl*r*(V^2)); % wing area
b = sqrt(AR * S); % wing span (tip to tip)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
run('MAE_154A_weight.m');

% Plane specs
rho = 0.0023769;  % air density (slug/ft^3)
V = 50;        % stall velocity (ft/s)
CL_MAX = 1.5;      % Cl max 
%W = 19;        % total weight of aircraft + avioinics/payload/etc. (lbs)
W = W_pick;

% Plane Parameter
%L_fuse = 5; % length of fuselage (ft) (defined in weight code)

% Wing Parameters
Arw = 7.4;         % aspect ratio
Sw = (W*2) / (CL_MAX*rho*(V^2)); % wing area
bw = sqrt(Arw * Sw); % wing span (tip to tip)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
%c = S/b;          % chord length 
c = 0.8; 

% Horizontal tail parameters

AR_ht = 4;        % horizontal tail aspect ratio (GUESSED!!)
V_ht = 0.5;       % horizontal tail volume coefficient (SEARCHED UP!!)
Lh = 4*c;         % distance from AC of wing to AC of horizontal tail
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
S_ht = (V_ht*S*c) / Lh; % Horizonatil tail area
b_ht = sqrt(AR_ht*S_ht);% Horizontal tail wing span
c_ht = S_ht / b_ht;     % Horizontal tail chord length

% Vertical tail parameters

AR_vt = 1.5;      % vertical tail aspect ratio (GUESSED!!)
V_vt = 0.05;      % volume coefficient vertical tail (SEARCHED UP!!)
Lv = 4*c;         % distance from AC of wing to AC of verticla tail
S_vt = (V_vt*S*c*b) / Lv; % Vertical tail area
b_vt = sqrt(AR_vt*S_vt);  % Vertcial tail span 
c_vt = S_vt / b_vt;       % Vertcial tail chord length 


% Structural Weights 

Wstruct = 4.5; % total weight of empty aircrafts (lbs)
Wf = 0.7;      % weight of fuselage (lbs)
Ww = 1.4;      % Weight of wing
Wht = 0.5233;  % weight of horizontal tail
Wvt = 0.3728;  % weight of vertical tail
Wlg = 1.5;     % weight of landing gear
Wpr = 0.5;     % weight of propellar

% locations of LE, AC, and CG of wing, htail, vtail, landing gear, & prop

x_wle = 0.75; % wing leading edge wrt nose (ft)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
S_th = (V_ht*Sw*c) / Lh; % Horizonatil tail area
Sth = 0.2 * Sw;
bth = sqrt(AR_ht*Sth);% Horizontal tail wing span
%c_ht = Sth / bth;     % Horizontal tail chord length
c_ht = 0.38;

% Vertical tail parameters ( the rest calculated further down)

AR_vt = 1.8;      % vertical tail aspect ratio (CONSTANT!!)
V_vt = 0.04;      % volume coefficient vertical tail (CONSTANT!!)
% Temporary!!
Lv = 4*c;         % distance from AC of wing to AC of verticla tail

% Structural Weights (use from other code)

% Wstruct = 4.3720; % total weight of empty aircrafts (lbs)
% Wf = 0.6465;      % weight of fuselage (lbs)
% Ww = 1.2398;      % Weight of wing
% Wht = 0.6363;  % weight of horizontal tail
% Wvt = 0.3493;  % weight of vertical tail
Wlg = 1.5;     % weight of landing gear
Wprop = .24;     % weight of propellar

% locations of LE, AC, and CG of wing, htail, vtail, landing gear, & prop

x_wle = 1; % wing leading edge wrt nose (ft)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
x_wac = x_wle + 0.25*c; % location of wing aerodynamic center wrt nose
x_wcg = x_wle + 0.35*c; % location of wing center of gravity wrt nose

x_hac = x_wac + Lh;        % location of horizontal tail AC wrt nose
x_hle = x_hac - 0.25*c_ht; % location of horizontal tail leading edge wrt nose
x_hcg = x_hle + 0.35*c_ht; % location of horizontal tail cg wrt nose

x_vac = x_wac + Lv;         % location of vertical tail AC wrt nose
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
x_vle = x_vac - 0.25*c_vt;  % location of vertical tail leading edge wrt nose
x_vcg = x_vle + 0.35*c_vt; % location of vertical tail cg wrt nose

x_fcg = 0.4*Lf; % CG of fuselage (estimated)
x_lgcg = 1.5; % CG of landing gear (estimated)
x_prcg = 0; % CG of propellar

% Total Center of gravity of structure (empty plane) !!!
W = [Ww, Wht, Wvt, Wf, Wlg, Wpr];
x = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_prcg];

x_cg_str = sum(W.*x) / sum(W) % Total structural center of gravity 

% Weight of Avionics (lbs)
Wcam = 1.87;  % weight of camera 
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
                            % the rest calculated further down

x_fcg = 0.4*L_fuse; % CG of fuselage (estimated)
x_lgcg = 2; % CG of landing gear (estimated)
x_propcg = -0.2; % CG of propellar

% Vertical Tail Parameters!! (contd)
W_novt = [Ww, Wht, Wf, Wlg, Wprop]; % weight without vertical tail
x_novt = [x_wcg, x_hcg, x_fcg, x_lgcg, x_propcg]; 

x_cg_noVT = sum(W_novt.*x_novt) / sum(W_novt); % CG withough vertical tail

lv = x_vac - x_cg_noVT   % distance from cg of plane to ac of v tail
Stv = V_vt*((Sw*bw)/lv)  % VERTICAL TAIL AREA
btv = sqrt(AR_vt*Stv)    % VERTICAl TAIL WINGSPAN
cv = Stv / btv           % VERTICAL TAIL CHORD (airfoil NACA 0012)

x_vle = x_vac - 0.25*cv;  % location of vertical tail leading edge wrt nose
x_vcg = x_vle + 0.35*cv; % location of vertical tail cg

% Total Center of gravity of structure (empty plane) !!!
W = [Ww, Wht, Wvt, Wf, Wlg, Wprop];
x = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_propcg];

x_cg_str = sum(W.*x) / sum(W); % Total structural center of gravity 

% Weight of Avionics (lbs)
Wcam = 0.95;  % weight of camera 
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
Wcomp = 0.10; % weight of computer
Wgps = 0.07;  % weight of GPS
Wbat = 0.44;  % weight of battery
Wserv = 0.06; % weight of servos
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
Weng = 3;     % weight of engine
=======
Weng = 1.76;     % weight of engine
Wfs = 0.25;   % weight of fuel system
>>>>>>> Stashed changes
=======
Weng = 1.76;     % weight of engine
Wfs = 0.25;   % weight of fuel system
>>>>>>> Stashed changes
=======
Weng = 1.76;     % weight of engine
Wfs = 0.25;   % weight of fuel system
>>>>>>> Stashed changes

% Weight of payload/fuel
Wbal = 0.66;   % weight of 1 waterballoon
Wpay = 3*Wbal; % weight of payload
%Wpay = 0;
Wfuel = 1.5;   % weight of initial fuel
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

% Locations of avionics wrt to nose
x_cam = 2;    % location of camera
x_comp = 3.5; % location of computer
x_gps = 3.5;  % location of GPS
x_bat = 2.5;  % location of battery
x_serv = 1.8; % location of servos
x_eng = 0.25; % location of engine

% Locations of payload/fuel wrt nose
x_1bal = 1.8;
X_2bal = 1.8;
x_3bal = 1.8;
x_pay = 1.8;  % location of payload
x_fuel = 1.8; % loaction of fuel 

% Total CG without payload or fuel!!
Wtot = [Ww, Wht, Wvt, Wf, Wlg, Wpr, Wcam, Wcomp, Wgps, Wbat, Wserv, Weng]; %weight of structure + avionics
x_tot  = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_prcg, x_cam, x_comp, x_gps, x_bat, x_serv, x_eng];

x_cg_tot = sum(Wtot.*x_tot)/sum(Wtot) 

% Initial CG with payload and fuel!!
W0 = [Ww, Wht, Wvt, Wf, Wlg, Wpr, Wcam, Wcomp, Wgps, Wbat, Wserv, Weng, Wpay, Wfuel]; % initial weight w/ payload and fuel
x_0  = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_prcg, x_cam, x_comp, x_gps, x_bat, x_serv, x_eng, x_pay, x_fuel];

x_cg_0 = sum(W0.*x_0)/sum(W0) % initial cg w/ payload and fuel

=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
Wfuel_half = 0.75;
Wfuel_quart = 0.375;

% Locations of avionics wrt to nose
x_cam = 0.78;    % location of camera
x_comp = 2.5; % location of computer
x_gps = 2.5;  % location of GPS
x_bat = 1;  % location of battery
x_serv = 1.25; % location of servos
x_eng = 0.25; % location of engine
x_fs = 1.25;   % location of fuel system 

% Locations of payload/fuel wrt nose
x_1bal = 1.25;
X_2bal = 1.25;
x_3bal = 1.25;
x_pay = 1.25;  % location of payload
x_fuel = 1.25; % loaction of fuel 

% Total CG without payload or fuel!!
Wdry = [W, Wcam, Wcomp, Wgps, Wbat, Wserv, Weng, Wfs]; %weight of structure + avionics
x_dry  = [x, x_cam, x_comp, x_gps, x_bat, x_serv, x_eng, x_fs];

x_cg_dry = sum(Wdry.*x_dry)/sum(Wdry)

% Initial CG with payload and fuel!!
W0 = [Wdry, Wpay, Wfuel]; % initial weight w/ payload and fuel
x_0  = [x_dry, x_pay, x_fuel];

x_cg_0 = sum(W0.*x_0)/sum(W0) % initial cg w/ payload and fuel

% CG with 1/2 fuel and all payload
W_fhalf_pay = [Wdry, Wpay, Wfuel_half]; % weight w/ payload and 1/2 fuel
x_fhalf_pay  = [x_dry, x_pay, x_fuel];

x_cg_fhalf_pay = sum(W_fhalf_pay.*x_fhalf_pay)/sum(W_fhalf_pay) % cg w/ payload and 1/2 fuel

% CG with 1/2 fuel and no payload
W_fhalf_nopay = [Wdry, Wfuel_half]; % weight w/ no payload and 1/2 fuel
x_fhalf_nopay  = [x_dry, x_fuel];

x_cg_fhalf_nopay = sum(W_fhalf_nopay.*x_fhalf_nopay)/sum(W_fhalf_nopay) % cg w/ no payload and 1/2 fuel

% CG with 1/4 fuel and no payload
W_fquart_nopay = [Wdry, Wfuel_quart]; % weight w/ no payload and 1/4 fuel
x_fquart_nopay  = [x_dry, x_fuel];

x_cg_fquart_nopay = sum(W_fquart_nopay.*x_fquart_nopay)/sum(W_fquart_nopay) % cg w/ no payload and 1/4 fuel

% CG with fuel and no payload
W_f_nopay = [Wdry, Wfuel]; % weight w/ no payload and all fuel
x_f_nopay  = [x_dry, x_fuel];

x_cg_f_nopay = sum(W_f_nopay.*x_f_nopay)/sum(W_f_nopay) % cg w/ no payload and all fuel

<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

