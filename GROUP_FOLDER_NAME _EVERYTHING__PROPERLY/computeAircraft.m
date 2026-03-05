function outputs = computeAircraft(inputs)

%% =====================================
% UNPACK INPUTS 
% =====================================
%fix
% --- BASE AIRCRAFT GEOMETRY & AERO ---
rho      = inputs.rho;
mu       = inputs.mu;
a        = inputs.a;
W        = inputs.W;
Arw      = inputs.Arw;
bw       = inputs.bw;
lambdaw  = inputs.lambdaw;

lambdath = inputs.lambdath;
Art      = inputs.Art;
e        = inputs.e;
cla      = inputs.cla;
downwash = inputs.downwash;
it       = inputs.it;
Cmacw    = inputs.Cmacw;
tau      = inputs.tau;
CL_MAX   = inputs.CL_MAX; 
Power    = inputs.Power;
EF       = inputs.EF;
L_fuse   = inputs.L_fuse;
W_fuse   = inputs.W_fuse;
H_fuse   = inputs.H_fuse;
D_fuse   = inputs.D_fuse;

% REMOVED FROM BASE: 
% x_wle    = inputs.x_wle;      % -> Moved to Adela's section (updated value)
% x_cg_dry = inputs.x_cg_dry;   % -> Adela's code calculates this dynamically now
% x_cg_half= inputs.x_cg_half;  % -> Adela's code calculates this dynamically now
% x_cg_0   = inputs.x_cg_0;     % -> Adela's code calculates this dynamically now


% --- CHANGES BY JONATHAN ---
c_vt     = inputs.c_vt; 
tc_v     = inputs.tc_v;     
xcm_v    = inputs.xcm_v; 
Lam_v    = inputs.Lam_v; 
V        = inputs.V; 

% REMOVED FROM JONATHAN:
% AR_vt  = inputs.AR_vt;        % -> Adela updated this to 1.5, kept in her section

% --- ADELA's STUFF ---

% Horizonatl Tail 
c_ht     = inputs.c_ht;
hac_ht   = inputs.hac_ht;
V_ht     = inputs.V_ht;

% Vertical Tail
hacv     = inputs.hacv;
AR_vt    = inputs.AR_vt;
V_vt     = inputs.V_vt;
Sv       = inputs.Sv;
bv       = inputs.bv;

% Wing
c          = inputs.c;
hac        = inputs.hac;

% Flight / Design Inputs for weight
Wguess     = inputs.Wguess;     
V_stall    = inputs.V_stall; 
N_load     = inputs.N_load;
tc         = inputs.tc;
V_max_kts  = inputs.V_max_kts;           

% Specific Component Weights
Wlg        = inputs.Wlg;
Wprop      = inputs.Wprop;
Weng       = inputs.Weng;
Neng       = inputs.Neng;
Wfuel      = inputs.Wfuel;
Wfs        = inputs.Wfs;
Wau        = inputs.Wau;
Wbal       = inputs.Wbal;
N_bal      = inputs.N_bal;
Wpay       = inputs.Wbal * inputs.N_bal;
Wcam       = inputs.Wcam;
Wcomp      = inputs.Wcomp;
Wgps       = inputs.Wgps;
Wbat       = inputs.Wbat;
Wserv      = inputs.Wserv;

% Component CG Locations
x_wle      = inputs.x_wle;
x_lgcg     = inputs.x_lgcg;
x_propcg   = inputs.x_propcg;
x_cam      = inputs.x_cam;
x_comp     = inputs.x_comp;
x_gps      = inputs.x_gps;
x_bat      = inputs.x_bat;
x_serv     = inputs.x_serv;
x_eng      = inputs.x_eng;
x_fs       = inputs.x_fs;
x_pay      = inputs.x_pay;
x_fuel     = inputs.x_fuel;
c_p_hp_hr  = inputs.c_p_hp_hr;

Sw = bw^2/Arw;
%% horizontal tail stuff 
Sth= 0.20*Sw ; 
bth=sqrt(Sw * Art) ; 

% Display Options
makePlots   = inputs.makePlots;

%% ADELAS MATLAB FILES WEiGHT and CG

% =====================================
% weight
% =====================================

% Initialize guess
W_current = inputs.Wguess;
Wto = zeros(1, 5); % Locked to 5 iterations as requested

for i = 1:10
    Wto(i) = W;
% Wing Weight 
    % lambdaw
    % tc
    % V_max_kts
    Delta       = 0*pi/180;
    Ww = 96.948*((W*N_load/10^5)^0.65*(Arw/cos(Delta))^0.57*(Sw/100)^0.61*((1+lambdaw)/(2*tc))^0.36*(1+V_max_kts/500)^0.5)^0.993;

% Fuselage Weight
     % L_fuse
     % W_fuse
     % D_fuse
    Wf=200*((W*N_load/10^5)^0.286*(L_fuse/10)^0.857*((W_fuse+D_fuse)/10)*(V_max_kts/100)^0.338)^1.1;

% Horizontal Tail Weight
    lh=35 / 12 + (.5 - hac) * c - (.5 - hac_ht) * c_ht; %ft %Distance from Wing MACto Tail MAC
    thr=c_ht*.12*12; %inches %horizontal tail max root thickness (chord * thick/chord)
    Wht=127*((W*N_load/10^5)^0.87*(Sth/100)^1.2*(lh/10)^0.483*(bth/thr)^0.5)^0.458;

% Landing Gear Weight
    % Wlg = 1.5;

% Weight of propellar 
    % Wprop = 0.24;

% Vertical Tail Weight
    tvr=c_vt*.12*12; %in %Vertical Tail Max Root Thickness (chord * thick/chord * in/ft)
    Wvt= (2)* 98.5*((W*N_load/10^5)^0.87*( (.5)* Sv/100)^1.2*( (.5)* bv/tvr)^0.5)^0.458;

% TOTAL STRUCTURAL WEIGHT
    Wstruct=Ww+Wf+Wht+Wvt+Wlg;

% Total Propulsion Unit (minus Fuel system) Weight
    % Weng
    % Neng
    Wp=2.575*(Weng)^0.922*Neng; %this equation likely over-estimates propulsion unit weight for a small UAV

% Fuel Weight
    % Wfuel = 1.5; %(lbs)

% Fuel System Weight
    % Wfs = .25;

% Surface Controls Weight
    Wsc=1.066*W^0.626;

% Avionics Weight (camera, servos, computer, GPS & Battery)
    % Wau=1.62;

% Payload Weight
    % Wpay=2;

%% TOTAL WEIGHT
    Wto(i)=Wstruct+Wp+Wfs+Wsc+Wpay+Wfuel;
    W = Wto(i);

    Ww_arr(i) = Ww;
    Wf_arr(i) = Wf;
    Wht_arr(i) = Wht;
    Wvt_arr(i) = Wvt;
    Wstruct_arr(i) = Wstruct;
    Wp_arr(i) = Wp;
    Wau_arr(i) = Wau;
end

% W at elbow!!
index = 3;            % should choose x=4
W_total = Wto(index);



%Sw = (W_total*2) / (CL_MAX*rho*(V_stall^2)); % wing area
%bw = sqrt(Arw * Sw); % wing span (tip to tip)

% Outputs of component weights @ index!!
Ww = Ww_arr(index);
Wf = Wf_arr(index);
Wht = Wht_arr(index);
Wvt = Wvt_arr(index);
Wstruct = Wstruct_arr(index);
Wp = Wp_arr(index);
Wau = Wau_arr(index);

% figure; grid on; hold on
% plot([Wguess Wto], '.-m')
% plot(index+1, W_total, 'or', 'MarkerSize', 10, 'LineWidth', 2)  % +1 bc Wguess is the first point

% disp('Final calculated total weight (5th iteration):');
% disp(Wto(5));

% =====================================
% Center of Gravity + V-tail calculations
% =====================================

W = W_total;
L_ht = 4*c;         % distance from AC of wing to AC of horizontal tail
Lv = 4*c;         % distance from AC of wing to AC of verticla tail

% locations of LE, AC, and CG of wing, htail, vtail, landing gear, & prop
x_wac = x_wle + 0.25*c; 
x_wcg = x_wle + 0.35*c;
x_hac = x_wac + L_ht;        
x_hle = x_hac - 0.25*c_ht; 
x_hcg = x_hle + 0.35*c_ht;
x_vac = x_wac + Lv;
x_fcg = 0.4*L_fuse; % CG of fuselage (estimated)

%% Vertical Tail (contd)
W_novt = [Ww, Wht, Wf, Wlg, Wprop];  % weight without vertical tail
x_novt = [x_wcg, x_hcg, x_fcg, x_lgcg, x_propcg];  

x_cg_noVT = sum(W_novt.*x_novt) / sum(W_novt); % cg without vertical tail

% Vertical Tail Parameters
l_vt = x_vac - x_cg_noVT; % distance from cg wihtout vt to ac of vt 
S_vt = V_vt*((Sw*bw)/l_vt); % VERTICAL TAIL AREA
b_vt = sqrt(AR_vt*S_vt);    % VERTICAL TAIL span
c_vt = S_vt / b_vt;           % VERTICAL TAIL chord (NACA 0012)

x_vle = x_vac - 0.25*c_vt;  % location of vertical tail leading edge wrt nose
x_vcg = x_vle + 0.35*c_vt; % location of vertical tail cg

%% Total Center of gravity of structure (empty plane) !!!
W = [Ww, Wht, Wvt, Wf, Wlg, Wprop];
x = [x_wcg, x_hcg, x_vcg, x_fcg, x_lgcg, x_propcg];

x_cg_str = sum(W.*x) / sum(W); % Total structural center of gravity 

Wpay = N_bal*Wbal; % weight of payload
Wfuel_half = 0.5*Wfuel;
Wfuel_quart = 0.25*Wfuel;

%% Locations of payload/fuel wrt nose
x_1bal = 1.25;
X_2bal = 1.25;
x_3bal = 1.25;
% x_pay = 1.25;  % location of payload
% x_fuel = 1.25; % loaction of fuel 

%% Total CG without payload or fuel!!
Wdry = [W, Wcam, Wcomp, Wgps, Wbat, Wserv, Weng, Wfs]; %weight of structure + avionics
x_dry  = [x, x_cam, x_comp, x_gps, x_bat, x_serv, x_eng, x_fs];
x_cg_dry = sum(Wdry.*x_dry)/sum(Wdry);

%% Initial CG with payload and fuel!!
W0 = [Wdry, Wpay, Wfuel]; % initial weight w/ payload and fuel
x_0  = [x_dry, x_pay, x_fuel];
x_cg_0 = sum(W0.*x_0)/sum(W0);

%% CG with 1/2 fuel and all payload
W_fhalf_pay = [Wdry, Wpay, Wfuel_half]; % weight w/ payload and 1/2 fuel
x_fhalf_pay  = [x_dry, x_pay, x_fuel];
x_cg_fhalf_pay = sum(W_fhalf_pay.*x_fhalf_pay)/sum(W_fhalf_pay);

%% CG with 1/2 fuel and no payload
W_fhalf_nopay = [Wdry, Wfuel_half]; % weight w/ no payload and 1/2 fuel
x_fhalf_nopay  = [x_dry, x_fuel];
x_cg_fhalf_nopay = sum(W_fhalf_nopay.*x_fhalf_nopay)/sum(W_fhalf_nopay);

%% CG with 1/4 fuel and no payload
W_fquart_nopay = [Wdry, Wfuel_quart]; % weight w/ no payload and 1/4 fuel
x_fquart_nopay  = [x_dry, x_fuel];
x_cg_fquart_nopay = sum(W_fquart_nopay.*x_fquart_nopay)/sum(W_fquart_nopay);

%% CG with fuel and no payload
W_f_nopay = [Wdry, Wfuel]; % weight w/ no payload and all fuel
x_f_nopay  = [x_dry, x_fuel];
x_cg_f_nopay = sum(W_f_nopay.*x_f_nopay)/sum(W_f_nopay);


%% (10) DISPLAY OPTIONS (directly from inputs)
makePlots            = inputs.makePlots;           
makeTable            = inputs.makeTable;           
makePrint_stability  = inputs.makePrint_stability; 
makePrint_tail_Volume = inputs.makePrint_tail_Volume;
%% Weight Adelas stuff 
% Map Adela's newly calculated half-fuel CG to the old variable name 
% so the stability math below doesn't break!
x_cg_half = x_cg_fhalf_pay; 

% (And if your code also looks for the old dry or initial CGs, map them too:)
x_cg_dry = x_cg_f_nopay; % Or whichever dry scenario you need
x_cg_0   = x_cg_0;






%% Fuselage
Swet_f = 2*(L_fuse*W_fuse + L_fuse*H_fuse + W_fuse*H_fuse);
d_fuse = W_fuse;  % for fuselage form factor

%% ========================================================================
% 2) PLANFORM GEOMETRY (ROOT/TIP/MAC)
% ========================================================================
Crw = (2*Sw)/(bw*(1 + lambdaw));
Ctw = lambdaw * Crw;
Cw  = (2/3) * Crw * ((1 + lambdaw + lambdaw^2)/(1 + lambdaw)); % wing MAC

Crth = (2*Sth)/(bth*(1 + lambdath));
Ctth = lambdath * Crth;
Cth  = (2/3) * Crth * ((1 + lambdath + lambdath^2)/(1 + lambdath)); % tail MAC
%% ========================================================================
% 3) 3D LIFT CURVE SLOPES
aw = cla / (1 + (cla/(pi*e*Arw)));
at = cla / (1 + (cla/(pi*e*Art)));

%% ========================================================================
% 4) AC LOCATIONS FROM GEOMETRY (NO GUESSING)
% Tail root TE flush with fuselage end; tail AC at 0.25*c behind tail LE.
% ========================================================================
hacw = 0.25;   % wing AC at 0.25 MAC
hact = 0.25;   % tail AC at 0.25 chord (root chord ref)

x_wac = x_wle + hacw*Cw;                 % wing AC from nose [ft]
x_tac = L_fuse - (1 - hact)*Crth;        % tail AC from nose [ft]
Lh    = x_tac - x_wac;                   % wing AC -> tail AC [ft]

x_ac_t_LE = x_tac - x_wle;               % tail AC measured from wing LE [ft]
h_ac_t    = x_ac_t_LE / Cw;              % nondim tail AC from wing LE [-]
%% ========================================================================
% 5) NEUTRAL POINT + STATIC MARGIN
% ========================================================================
hn = (hacw + h_ac_t*((Sth/Sw)*(at/aw))*(1-downwash)) / (1 + (Sth/Sw)*(at/aw)*(1-downwash));

hcg_dry  = (x_cg_dry  - x_wle)/Cw;
hcg_half = (x_cg_half - x_wle)/Cw;
hcg_0    = (x_cg_0    - x_wle)/Cw;

SM_dry  = hn - hcg_dry;
SM_half = hn - hcg_half;
SM_0    = hn - hcg_0;
%% ========================================================================
% 6) TAIL ARM (CG -> TAIL AC) AND TAIL VOLUME VH
% ========================================================================
lt_dry  = x_tac - x_cg_dry;
lt_half = x_tac - x_cg_half;
lt_0    = x_tac - x_cg_0;

VH_dry  = (Sth*lt_dry )/(Sw*Cw);
VH_half = (Sth*lt_half)/(Sw*Cw);
VH_0    = (Sth*lt_0   )/(Sw*Cw);

%% ========================================================================
% 7) LONGITUDINAL COEFFICIENTS (KEEP ALL)
% Use “half fuel + payload” as your design case (you can swap later).
% ========================================================================
hcg = hcg_half;
VH  = VH_half;

% Aircraft lift slope
CLalpha = aw + at*(Sth/Sw)*(1 - downwash);

% --- These three are the ones you asked to KEEP explicit ---
CL0 = -at*(Sth/Sw)*it;      % lift at alpha=0 due to incidence
CLi = -at*(Sth/Sw);         % dCL/d(it)
CMi =  at*VH;               % dCm/d(it)

% Pitching moment slope (static stability)
Cmalpha = -CLalpha*(hn - hcg);

% Zero-alpha pitching moment about CG
Cm0 = Cmacw + VH*at*it;

% Elevator derivatives
CLdeltae = at*(Sth/Sw)*tau;
Cmdeltae = -at*VH*tau;

%% ========================================================================
% 8) VELOCITY SWEEP
% ========================================================================
q = 0.5*rho.*V.^2;
M = V./a;

% Required CL for level flight (L=W)
CL = W_total ./ (q.*Sw);

% Solve for alpha, alpha_t (first-pass linear)
alpha   = CL ./ aw;                         % [rad] wing-only first-pass
alpha_t = (1 - downwash).*alpha + it;       % [rad] tail AoA model

% Wing lift coefficient contribution
CLw = aw .* alpha;

% Elevator deflection needed for trim at each V (kept explicit for team)
% deltae(CL) = -((Cm0*CLalpha) + (Cmalpha*CL)) / ((CLalpha*Cmdeltae) - (Cmalpha*CLdeltae))
deltae = -((Cm0*CLalpha) + (Cmalpha.*CL)) ./ ((CLalpha*Cmdeltae) - (Cmalpha*CLdeltae)); % [rad]

% Tail lift coefficient using YOUR professor-style form:
% CLt = CL0 + at*(alpha_t - it) + CLdeltae*deltae
CLt = CL0 + at*(alpha_t - it) + CLdeltae.*deltae;
%% ========================================================================
% 9) PARASITE DRAG (Raymer-style build-up) + INCLUDE VERTICAL TAIL
% ========================================================================
% Wetted areas (first-pass)
Swet_w = 2*Sw;
Swet_t = 2*Sth;
Swet_v = 2*S_vt;

% Reynolds numbers (needed internally for Cf; we just won’t print them)
Re_w = rho.*V.*Cw    ./ mu;
Re_t = rho.*V.*Cth   ./ mu;
Re_v = rho.*V.*c_vt    ./ mu;
Re_f = rho.*V.*L_fuse./ mu;

% Skin friction (Raymer-style)
Cf_w = 0.455 ./ ((log10(Re_w)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_t = 0.455 ./ ((log10(Re_t)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_v = 0.455 ./ ((log10(Re_v)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_f = 0.455 ./ ((log10(Re_f)).^2.58 .* (1 + 0.144*M.^2).^0.65);

% Surface form factor (same expression you’ve been using)
tc  = 0.12;
xcm = 0.30;
Lam = 0;
K_surf = (1 + (0.6/xcm)*tc + 100*tc^4) .* (1.34*M.^0.18 .* (cos(Lam)).^0.28);

% Vertical tail form factor uses its own tc/xcm/Lam (kept separate)
K_v = (1 + (0.6/xcm_v)*tc_v + 100*tc_v^4) .* (1.34*M.^0.18 .* (cos(Lam_v)).^0.28);

% Interference factors
Q_w = 1.0;
Q_t = 1.05;
Q_v = 1.05;
Q_f = 1.0;

% Fuselage form factor
f   = L_fuse/d_fuse;
K_f = 0.9 + 5/f^1.5 + f/400;

% Misc allowances
CD_misc = 0.002;
CD_LP   = 0.002;

% Total parasite CD
CDp = ...
    K_surf.*Q_w.*Cf_w.*(Swet_w/Sw) + ...
    K_surf.*Q_t.*Cf_t.*(Swet_t/Sw) + ...
    K_v   .*Q_v.*Cf_v.*(Swet_v/Sw) + ...
    K_f   .*Q_f.*Cf_f.*(Swet_f/Sw) + ...
    CD_misc + CD_LP;
%% ========================================================================
% 10) INDUCED DRAG + TOTAL DRAG + L/D
% ========================================================================
CDi_w = (CLw.^2) ./ (pi*Arw*e);
CDi_t = (Sth/Sw) .* (CLt.^2) ./ (pi*Art*e);
CDi_total = CDi_w + CDi_t;

CD_total = CDp + CDi_total;

Dp = q .* Sw .* CDp;
Di = q .* Sw .* CDi_total;
D  = q .* Sw .* CD_total;

L  = q .* Sw .* CL;     % should match W
LD = L ./ D;

% Power required (hp)
Power_required_hp = (D.*V) / 550;
Power_available_hp = Power * ones(size(V))*EF;
%% ========================================================================
% 11) FIND (L/D)max AND PRINT KEY VALUES
% ========================================================================
[LDmax, idx] = max(LD);

V_LDmax        = V(idx);
CL_LDmax       = CL(idx);
alpha_LDmax_deg   = rad2deg(alpha(idx));
alpha_t_LDmax_deg = rad2deg(alpha_t(idx));
deltae_LDmax_deg  = rad2deg(deltae(idx));



%% Stall speed (velocity where CL = CL_MAX)
V_stall = sqrt(2 * W_total / (rho * Sw * CL_MAX));

%% Maximum speed (velocity where Power Required <= Power Available)
V_max_idx = find(Power_required_hp <= Power_available_hp, 1, 'last');
if ~isempty(V_max_idx)
    V_max = V(V_max_idx);
else
    V_max = NaN;
end

%% Rate of climb (ft/s)
ROC = (Power_available_hp - Power_required_hp) * 550 ./ W_total;  % ft*lbf/s / lb = ft/s

% Rate of climb at stall
if ~isnan(V_stall)
    ROC_stall = ROC(find(V >= V_stall, 1, 'first'));
else
    ROC_stall = NaN;
end

%% 
% Inputs:
% eta_pr = propeller efficiency
% c_p    = specific fuel consumption (1/hr or 1/s depending on units)
% CL     = lift coefficient
% CD     = drag coefficient
% rho    = air density (slug/ft^3 for Imperial)
% S      = wing area (ft^2)
% Wi     = initial weight (lbf)
% Wf     = final weight (lbf)
%% ========================================================================
% 12) MAXIMUM ENDURANCE (Breguet Equation for Propeller Aircraft)
% ========================================================================
% Max endurance occurs where (CL^1.5 / CD) is maximized.
power_index = (CL.^1.5) ./ CD_total; 
[~, end_idx] = max(power_index);

% Grab the CL and CD at that specific flight condition
CL_endurance = CL(end_idx);
CD_endurance = CD_total(end_idx);

% Convert c_p from [lb/(hp*hr)] to standard base units [1/ft]
% 1 hp = 550 ft-lbf/s, 1 hour = 3600 seconds
c_p_base = c_p_hp_hr / (550 * 3600);

% Define Initial (Wi) and Final (Wf) weights in lbs
Wi = W_total;
Wf = W_total - Wfuel;

% Calculate Endurance in seconds
Endurance_sec = (EF * (CL_endurance^1.5) / (c_p_base * CD_endurance)) ...
                * sqrt(2 * rho * Sw) * (1/sqrt(Wf) - 1/sqrt(Wi));

% Convert to minutes
Endurance_min = Endurance_sec / 60;

%% ========================================================================
% ) OUTPUT STORAGE
outputs = struct;
outputs.V = V;
outputs.CL = CL;
outputs.alpha = alpha;
outputs.alpha_t = alpha_t;
outputs.CLw = CLw;
outputs.CLt = CLt;
outputs.CDp = CDp;
outputs.CDi_w = CDi_w;
outputs.CDi_t = CDi_t;
outputs.CD_total = CD_total;
outputs.Dp = Dp;
outputs.Di = Di;
outputs.D = D;
outputs.LD = LD;
outputs.V_LDmax = V_LDmax;
outputs.LDmax = LDmax;
outputs.Power_required_hp = Power_required_hp;
outputs.Power_available_hp = Power_available_hp;
outputs.Re_w = Re_w;
outputs.Re_t = Re_t;
outputs.Re_f = Re_f;
outputs.lt_dry = lt_dry;
outputs.lt_half = lt_half;
outputs.lt_0 = lt_0;
outputs.VH_dry = VH_dry;
outputs.VH_half = VH_half;
outputs.VH_0 = VH_0;
outputs.hn = hn;
outputs.hcg_dry = hcg_dry;
outputs.hcg_half = hcg_half;
outputs.hcg_0 = hcg_0;
outputs.SM_dry = SM_dry;    % static margin
outputs.SM_half = SM_half;  % static margin
outputs.SM_0 = SM_0;        % static margin 
outputs.V_stall = V_stall;
outputs.V_max   = V_max;
outputs.ROC     = ROC;
outputs.ROC_stall = ROC_stall;   
%% new from jonathan
%
outputs.deltae = deltae ; 
outputs.CLalpha  = CLalpha;
outputs.CL0      = CL0;
outputs.CLi      = CLi;
outputs.CMi      = CMi;
outputs.Cmalpha  = Cmalpha;
outputs.Cm0      = Cm0;
outputs.CLdeltae = CLdeltae;
outputs.Cmdeltae = Cmdeltae;
outputs.alpha_deg = rad2deg(alpha);
outputs.alpha_t_deg = rad2deg(alpha_t);
outputs.deltae_deg = rad2deg(deltae);

outputs.V_LDmax            = V_LDmax;
outputs.LDmax              = LDmax;
outputs.CL_LDmax           = CL_LDmax;
outputs.alpha_LDmax_deg    = alpha_LDmax_deg;
outputs.alpha_t_LDmax_deg  = alpha_t_LDmax_deg;
outputs.deltae_LDmax_deg   = deltae_LDmax_deg;

%% Adela outputs 
%outputs.W_total   = Wto(index);%what is this      
outputs.W_total   = W_total;   
outputs.W_struct  = Wstruct;     
outputs.W_wing    = Ww;          
outputs.W_fuse    = Wf;          
outputs.W_htail   = Wht;         
outputs.W_vtail   = Wvt;         

% Dynamically Sized Geometry
outputs.S_wing    = Sw;           
outputs.b_wing    = bw;           
outputs.S_htail   = Sth;        
outputs.S_vtail   = S_vt;        

% Center of Gravity Locations [ft from nose]
outputs.x_cg_str           = x_cg_str;          % CG of empty structure
outputs.x_cg_dry           = x_cg_dry;          % CG of structure + avionics (Zero Fuel)
outputs.x_cg_0             = x_cg_0;            % Full payload, full fuel
outputs.x_cg_fhalf_pay     = x_cg_fhalf_pay;    % Full payload, 1/2 fuel
outputs.x_cg_fhalf_nopay   = x_cg_fhalf_nopay;  % No payload, 1/2 fuel
outputs.x_cg_fquart_nopay  = x_cg_fquart_nopay; % No payload, 1/4 fuel
outputs.x_cg_f_nopay       = x_cg_f_nopay;      % No payload, full fuel





%%
outputs.ROC_stall = ROC_stall;
outputs.Endurance_min = Endurance_min;

fprintf('Total Weight: %.2f lb\n', outputs.W_total);

if makePlots
    % Power
    figure; hold on
    plot(V, Power_Required_hp,'b-','LineWidth',2)
    plot(V, P_available_hp,'r--','LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('Power [hp]')
    title('Power Required vs Velocity'); legend('Power Req','Power Avail'); grid on

    % Drag
    figure; hold on
    plot(V,Dp,'LineWidth',2); plot(V,Di,'LineWidth',2); plot(V,D,'LineWidth',2)
    plot(V_LDmax,D(idx),'ko','MarkerSize',8,'LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('Drag [lb]')
    title('Drag Forces vs Velocity'); legend('Dp','Di','D','(L/D)_{max}'); grid on

    % L/D
    figure; hold on
    plot(V,LD,'LineWidth',2); plot(V_LDmax,LDmax,'ko','MarkerSize',8,'LineWidth',2)
    xlabel('Velocity [ft/s]'); ylabel('L/D'); title('L/D vs Velocity'); legend('L/D','(L/D)_{max}'); grid on
end

if makeTable
    alpha_deg = rad2deg(alpha); alpha_t_deg = rad2deg(alpha_t);
    T = table(V(:), CL(:), alpha_deg(:), alpha_t_deg(:), CLw(:), CLt(:), ...
        CDp(:), CDi_w(:), CDi_t(:), CD_total(:), Dp(:), Di(:), D(:), LD(:), ...
        Re_w(:), Re_t(:), Re_f(:), ...
        'VariableNames', {'V_ft_s','CL','alpha_deg','alpha_t_deg','CLw','CLt','CDp','CDi_w','CDi_t','CD_total','Dp_lb','Di_lb','D_lb','L_over_D','Re_w','Re_t','Re_f'});
    disp(T)
end

if makePrint_tail_Volume
    fprintf("\n==================== TAIL VOLUME ====================\n");
    fprintf("Tail AC station x_tac = %.4f ft\n", x_tac);
    fprintf("lt_str = %.4f ft --> VH_str = %.4f\n", lt_str, VH_str);
    fprintf("lt_tot = %.4f ft --> VH_tot = %.4f\n", lt_tot, VH_tot);
    fprintf("lt_0   = %.4f ft --> VH_0   = %.4f\n", lt_0,   VH_0);
end    

if makePrint_stability
    fprintf("\n==================== STABILITY ====================\n");
    fprintf("Neutral point hn = %.4f (fraction of MAC from wing LE)\n", hn);
    fprintf("hcg_str = %.4f  --> SM_str = %.4f (%.1f%% MAC)\n", hcg_str, SM_str, 100*SM_str);
    fprintf("hcg_tot = %.4f  --> SM_tot = %.4f (%.1f%% MAC)\n", hcg_tot, SM_tot, 100*SM_tot);
    fprintf("hcg_0   = %.4f  --> SM_0   = %.4f (%.1f%% MAC)\n", hcg_0, SM_0, 100*SM_0);
end
