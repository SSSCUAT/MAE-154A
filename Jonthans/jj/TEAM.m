%% ========================================================================
% MAE 154A - ONE SCRIPT (FIRST-PASS)
% Longitudinal stability + trim coefficients + drag build-up + V-sweep
% Keeps ALL coefficients (CL0, CLi, CMi, etc.) for team clarity.
% ========================================================================

clear; clc; close all;

% ========================================================================
% 0) ATMOSPHERE (SEA LEVEL)
% ========================================================================
rho = 0.0023769;        % [slug/ft^3]
mu  = 3.737e-7;         % [slug/(ft*s)]
a   = 1116.45;          % [ft/s]

% ========================================================================
% 1) AIRCRAFT INPUTS (EDIT WHEN DESIGN CHANGES)
% ========================================================================

% --- Weight / power ---
W_total     = 17.2;             % [lb]
Power = 2.8 ;            % [hp] power available (flat line first-pass)

% --- Wing geometry ---
Arw     = 6.5038;          % [-]
bw      = 4.5634;        % [ft]
Sw      = bw^2/Arw;  % [ft^2]
lambdaw = 0.40;         % [-]
display(Sw) 
% --- Horizontal tail geometry ---
Sth      = Sw*0.20;     % [ft^2]
Art      = 3.1945;         % [-]
bth      = sqrt(Sth*Art); % [ft]
lambdath = 0.50;        % [-]
display(bth)
% --- Vertical tail geometry (DEFINE THESE) ---
% NOTE: Put real values here when you have them.
S_vt   = 0.30;            % [ft^2]  
AR_vt  = 1.8;             % [-]  
b_vt = 0.4077735426 ; 
c_vt  = 0.35;            % [ft]    
tc_v = 0.12;            % [-]     
xcm_v = 0.30;           % [-]     
Lam_v = 0;              % [rad]   

% --- Aero assumptions ---
e        = 0.85;        % [-]
cla      = 5.73;        % [1/rad]
downwash = 0.25;        % [-] dε/dα
tau      = 0.5;         % [-] elevator effectiveness
Cmacw    = -0.05;       % [-] wing Cm about wing AC 

% --- FIXED TAIL INCIDENCE (FIRST PASS) ---
it = -1*pi/180;          % [rad] 

% --- Locations (from nose) ---
x_wle     = 1.0786;        % [ft] wing leading edge location
x_cg_dry  = 1.3843;     % [ft]
x_cg_half = 1.3492;     % [ft] 
x_cg_0    = 1.3425;     % [ft]

% .Fuselage 
L_fuse = 3.5270;           % [ft]
W_fuse = 8/12;          % [ft]
H_fuse = 8/12;          % [ft]
Swet_f = 2*(L_fuse*W_fuse + L_fuse*H_fuse + W_fuse*H_fuse); % [ft^2]
d_fuse = W_fuse;        % [ft] equiv diameter for fineness ratio

% ========================================================================
% 2) PLANFORM GEOMETRY (ROOT/TIP/MAC)
% ========================================================================
Crw = (2*Sw)/(bw*(1 + lambdaw));
Ctw = lambdaw * Crw;
Cw  = (2/3) * Crw * ((1 + lambdaw + lambdaw^2)/(1 + lambdaw)); % wing MAC

Crth = (2*Sth)/(bth*(1 + lambdath));
Ctth = lambdath * Crth;
Cth  = (2/3) * Crth * ((1 + lambdath + lambdath^2)/(1 + lambdath)); % tail MAC

fprintf("\n==================== GEOMETRY ====================\n");
fprintf("Wing: Cr=%.4f ft, Ct=%.4f ft, MAC Cw=%.4f ft\n", Crw, Ctw, Cw);
fprintf("HT:   Cr=%.4f ft, Ct=%.4f ft, MAC Cth=%.4f ft\n", Crth, Ctth, Cth);
fprintf("Fuse: L=%.2f ft, W=%.3f ft, H=%.3f ft, Swet_f=%.3f ft^2\n", L_fuse, W_fuse, H_fuse, Swet_f);

% ========================================================================
% 3) 3D LIFT CURVE SLOPES
% ========================================================================
aw = cla / (1 + (cla/(pi*e*Arw)));   % wing 3D slope [1/rad]
at = cla / (1 + (cla/(pi*e*Art)));   % tail 3D slope [1/rad]

fprintf("\n==================== LIFT SLOPES ====================\n");
fprintf("aw = %.4f 1/rad\n", aw);
fprintf("at = %.4f 1/rad\n", at);

% ========================================================================
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

fprintf("\n==================== REFERENCE LOCATIONS ====================\n");
fprintf("x_wac = %.4f ft, x_tac = %.4f ft, Lh = %.4f ft\n", x_wac, x_tac, Lh);
fprintf("h_ac_t = %.4f (tail AC from wing LE / Cw)\n", h_ac_t);

% ========================================================================
% 5) NEUTRAL POINT + STATIC MARGIN
% ========================================================================
hn = (hacw + h_ac_t*((Sth/Sw)*(at/aw))*(1-downwash)) / (1 + (Sth/Sw)*(at/aw)*(1-downwash));

hcg_dry  = (x_cg_dry  - x_wle)/Cw;
hcg_half = (x_cg_half - x_wle)/Cw;
hcg_0    = (x_cg_0    - x_wle)/Cw;

SM_dry  = hn - hcg_dry;
SM_half = hn - hcg_half;
SM_0    = hn - hcg_0;

fprintf("\n==================== STABILITY ====================\n");
fprintf("hn = %.4f\n", hn);
fprintf("SM_dry  = %.6f (%.1f%% MAC)\n", SM_dry,  100*SM_dry);
fprintf("SM_half = %.6f (%.1f%% MAC)\n", SM_half, 100*SM_half);
fprintf("SM_0    = %.6f (%.1f%% MAC)\n", SM_0,    100*SM_0);
fprintf("hcg_dry = %.4f\n", hcg_dry);
fprintf("hcg_half = %.4f\n", hcg_half);
fprintf("hcg_0 = %.4f\n", hcg_0);

% ========================================================================
% 6) TAIL ARM (CG -> TAIL AC) AND TAIL VOLUME VH
% ========================================================================
lt_dry  = x_tac - x_cg_dry;
lt_half = x_tac - x_cg_half;
lt_0    = x_tac - x_cg_0;

VH_dry  = (Sth*lt_dry )/(Sw*Cw);
VH_half = (Sth*lt_half)/(Sw*Cw);
VH_0    = (Sth*lt_0   )/(Sw*Cw);

fprintf("\n==================== TAIL ARM & VOLUME ====================\n");
fprintf("lt_half = %.4f ft, VH_half = %.4f\n", lt_half, VH_half);

% ========================================================================
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

% ========================================================================
% 8) VELOCITY SWEEP
% ========================================================================
V = linspace(30, 300, 500);     % [ft/s]
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

% Tail lift coefficient 
% CLt = CL0 + at*(alpha_t - it) + CLdeltae*deltae
CLt = CL0 + at*(alpha_t - it) + CLdeltae.*deltae;

% ========================================================================
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

% ========================================================================
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
% Power_required_hp = (D.*V) / 550;
% Power_available_hp = Power * ones(size(V));


% 
% % Power required (hp)
Power_required_hp = (D.*V) / 550;
% 
% % ========================================================================
% % PROP EFFICIENCY MODEL (16x8 prop at 8000 RPM - first pass)
% % Speed data taken from your prop plot in mph, then converted to ft/s
% % ========================================================================
% V_mph_data = [0 10 20 30 40 50 60 70 80];
% eta_data   = [0.00 0.25 0.48 0.63 0.72 0.74 0.70 0.55 0.00];
% 
% % Convert mph to ft/s
% V_prop_data = V_mph_data * 1.46667;
% 
% % Interpolate prop efficiency onto aircraft velocity grid
% eta_p = interp1(V_prop_data, eta_data, V, 'pchip', 0);
% 
% % Use your already-defined available shaft power
% P_shaft_hp = Power;
% 
% % Power and thrust available
% Power_available_hp = eta_p .* P_shaft_hp;
% Thrust_available   = 550 .* Power_available_hp ./ V;







% ========================================================================
% MULTI-RPM POWER AVAILABLE MODEL FOR 16x8 PROP
% All speeds below entered in mph from your prop plots
% Then converted to ft/s for use with aircraft code
% ========================================================================

mph_to_fts = 1.46667;

% -----------------------------
% Common aircraft reference speeds
% -----------------------------
V_stall_ref      = 55;   % ft/s
V_endurance_ref  = 60;   % ft/s  (max CL^(3/2)/CD point)
% top speed achievable will be computed from max envelope later

% -----------------------------
% Shaft power assumption
% IMPORTANT:
% Use full engine power for max-performance comparisons
% -----------------------------
P_shaft_hp = 2.8;

% ========================================================================
% RPM = 2000
% ========================================================================
V2000_mph  = [0 2 4 6 8 10 12 14 16 18 20];
eta2000    = [0.00 0.12 0.22 0.31 0.39 0.47 0.55 0.60 0.60 0.56 0.14];

V2000 = V2000_mph * mph_to_fts;
eta_p_2000 = interp1(V2000, eta2000, V, 'pchip', 0);
Pavail_2000 = eta_p_2000 .* P_shaft_hp;

% ========================================================================
% RPM = 3000
% ========================================================================
V3000_mph  = [0 3 6 9 12 15 18 21 24 27 30];
eta3000    = [0.00 0.18 0.32 0.44 0.54 0.60 0.64 0.65 0.63 0.50 0.22];

V3000 = V3000_mph * mph_to_fts;
eta_p_3000 = interp1(V3000, eta3000, V, 'pchip', 0);
Pavail_3000 = eta_p_3000 .* P_shaft_hp;

% ========================================================================
% RPM = 4000
% ========================================================================
V4000_mph  = [0 5 10 15 20 25 30 35 38 40 41];
eta4000    = [0.00 0.18 0.38 0.52 0.62 0.67 0.67 0.61 0.45 0.28 0.00];

V4000 = V4000_mph * mph_to_fts;
eta_p_4000 = interp1(V4000, eta4000, V, 'pchip', 0);
Pavail_4000 = eta_p_4000 .* P_shaft_hp;

% ========================================================================
% RPM = 5000
% ========================================================================
V5000_mph  = [0 5 10 15 20 25 30 35 40 45 48 50 51];
eta5000    = [0.00 0.18 0.34 0.46 0.56 0.64 0.68 0.70 0.69 0.63 0.48 0.30 0.00];

V5000 = V5000_mph * mph_to_fts;
eta_p_5000 = interp1(V5000, eta5000, V, 'pchip', 0);
Pavail_5000 = eta_p_5000 .* P_shaft_hp;

% ========================================================================
% RPM = 6000
% ========================================================================
V6000_mph  = [0 5 10 15 20 25 30 35 40 45 50 55 58 60 61];
eta6000    = [0.00 0.18 0.30 0.42 0.52 0.60 0.65 0.69 0.71 0.71 0.68 0.60 0.35 0.12 0.00];

V6000 = V6000_mph * mph_to_fts;
eta_p_6000 = interp1(V6000, eta6000, V, 'pchip', 0);
Pavail_6000 = eta_p_6000 .* P_shaft_hp;

% ========================================================================
% RPM = 7000
% ========================================================================
V7000_mph  = [0 10 20 30 40 50 55 60 65 68 70 71];
eta7000    = [0.00 0.25 0.48 0.60 0.68 0.73 0.72 0.62 0.37 0.12 0.00 0.00];

V7000 = V7000_mph * mph_to_fts;
eta_p_7000 = interp1(V7000, eta7000, V, 'pchip', 0);
Pavail_7000 = eta_p_7000 .* P_shaft_hp;

% ========================================================================
% RPM = 8000
% ========================================================================
V8000_mph  = [0 10 20 30 40 50 60 70 75 78 80 82];
eta8000    = [0.00 0.25 0.40 0.55 0.64 0.71 0.74 0.72 0.64 0.40 0.10 0.00];

V8000 = V8000_mph * mph_to_fts;
eta_p_8000 = interp1(V8000, eta8000, V, 'pchip', 0);
Pavail_8000 = eta_p_8000 .* P_shaft_hp;

% ========================================================================
% RPM = 9000
% ========================================================================
V9000_mph  = [0 10 20 30 40 50 60 70 80 85 88 90 92];
eta9000    = [0.00 0.25 0.40 0.55 0.64 0.71 0.74 0.75 0.70 0.57 0.40 0.10 0.00];

V9000 = V9000_mph * mph_to_fts;
eta_p_9000 = interp1(V9000, eta9000, V, 'pchip', 0);
Pavail_9000 = eta_p_9000 .* P_shaft_hp;

% ========================================================================
% RPM = 10000
% ========================================================================
V10000_mph  = [0 10 20 30 40 50 60 70 80 90 95 100 103];
eta10000    = [0.00 0.19 0.36 0.49 0.59 0.67 0.72 0.75 0.75 0.68 0.59 0.42 0.00];

V10000 = V10000_mph * mph_to_fts;
eta_p_10000 = interp1(V10000, eta10000, V, 'pchip', 0);
Pavail_10000 = eta_p_10000 .* P_shaft_hp;

% ========================================================================
% POWER AVAILABLE ENVELOPE
% ========================================================================
Pavail_all = [Pavail_2000;
              Pavail_3000;
              Pavail_4000;
              Pavail_5000;
              Pavail_6000;
              Pavail_7000;
              Pavail_8000;
              Pavail_9000;
              Pavail_10000];

Pavail_env = max(Pavail_all, [], 1);

% Compute top speed achievable from envelope intersection
idx_top = find(Pavail_env >= Power_required_hp, 1, 'last');
if ~isempty(idx_top)
    V_top_achievable = V(idx_top);
else
    V_top_achievable = NaN;
end













% ========================================================================
% 10.5) ENDURANCE METRIC FOR PROPELLER AIRCRAFT
% Max endurance occurs at max(CL^(3/2)/CD)
% ========================================================================

% Endurance metric from Breguet prop endurance condition
CL32_over_CD = (CL.^(3/2)) ./ CD_total;

% Find maximum endurance point
[CL32CD_max, idx_endurance] = max(CL32_over_CD);

V_endurance        = V(idx_endurance);
CL_endurance       = CL(idx_endurance);
CD_endurance       = CD_total(idx_endurance);
alpha_endurance_deg   = rad2deg(alpha(idx_endurance));
alpha_t_endurance_deg = rad2deg(alpha_t(idx_endurance));
deltae_endurance_deg  = rad2deg(deltae(idx_endurance));
%Power_req_endurance   = Power_required_hp(idx_endurance);

fprintf("\n==================== ENDURANCE CONDITION ====================\n");
fprintf("Max (CL^(3/2)/CD) = %.4f at V = %.2f ft/s\n", CL32CD_max, V_endurance);
fprintf("CL at endurance condition = %.4f\n", CL_endurance);
fprintf("CD at endurance condition = %.4f\n", CD_endurance);
fprintf("alpha at endurance condition = %.2f deg\n", alpha_endurance_deg);
fprintf("alpha_t at endurance condition = %.2f deg\n", alpha_t_endurance_deg);
fprintf("deltae at endurance condition = %.2f deg\n", deltae_endurance_deg);
%fprintf("Power required at endurance condition = %.4f hp\n", Power_req_endurance);










% ========================================================================
% 11) FIND (L/D)max AND PRINT KEY VALUES
% ========================================================================
[LDmax, idx] = max(LD);

V_LDmax        = V(idx);
CL_LDmax       = CL(idx);
alpha_LDmax_deg   = rad2deg(alpha(idx));
alpha_t_LDmax_deg = rad2deg(alpha_t(idx));
deltae_LDmax_deg  = rad2deg(deltae(idx));


CL_max = 1.5;  
V_stall = sqrt( (2*W_total) / (rho*Sw*CL_max) );


fprintf("\n==================== PERFORMANCE ====================\n");
fprintf("(L/D)max = %.3f at V = %.2f ft/s\n", LDmax, V_LDmax);
fprintf("CL(L/Dmax) = %.4f\n", CL_LDmax);
fprintf("alpha(L/Dmax) = %.2f deg\n", alpha_LDmax_deg);
fprintf("alpha_t(L/Dmax) = %.2f deg\n", alpha_t_LDmax_deg);
fprintf("deltae(L/Dmax) = %.2f deg\n", deltae_LDmax_deg);
fprintf("V_stall = %.2f ft/s\n", V_stall);




% ========================================================================
% 11.5) DRAG DERIVATIVES AROUND TRIM CONDITION
% Compute CD_alpha and CD_deltae about the trim point at (L/D)max
% ========================================================================

% -----------------------------
% Trim point = point at (L/D)max
% -----------------------------
i_trim = idx;

alpha_trim   = alpha(i_trim);      % [rad]
deltae_trim  = deltae(i_trim);     % [rad]
alpha_t_trim = alpha_t(i_trim);    % [rad]

CD_trim   = CD_total(i_trim);
CDp_trim  = CDp(i_trim);
CDi_w_trim = CDi_w(i_trim);

fprintf("\n==================== DRAG DERIVATIVES AT TRIM ====================\n");
fprintf("Trim condition taken at (L/D)max:\n");
fprintf("V_trim      = %.2f ft/s\n", V(i_trim));
fprintf("alpha_trim  = %.4f rad  (%.2f deg)\n", alpha_trim, rad2deg(alpha_trim));
fprintf("deltae_trim = %.4f rad  (%.2f deg)\n", deltae_trim, rad2deg(deltae_trim));
fprintf("CD_trim     = %.6f\n", CD_trim);

%% ------------------------------------------------------------------------
% A) CD_alpha from local linear fit of CD vs alpha
% Professor said: take a small segment and fit a line
% -------------------------------------------------------------------------

% choose a small neighborhood around trim: +/- 1 deg
dalpha_fit = deg2rad(1.0);

fit_mask = (alpha >= alpha_trim - dalpha_fit) & (alpha <= alpha_trim + dalpha_fit);

% make sure enough points exist
if sum(fit_mask) >= 2
    p_CD_alpha = polyfit(alpha(fit_mask), CD_total(fit_mask), 1);
    CD_alpha = p_CD_alpha(1);    % [1/rad]
    
    fprintf("CD_alpha (local fit) = %.6f 1/rad\n", CD_alpha);
else
    CD_alpha = NaN;
    warning('Not enough points in alpha window to compute CD_alpha.');
end

%% ------------------------------------------------------------------------
% B) CD_deltae from small perturbation in elevator deflection at trim
% Keep trim flight condition fixed and perturb deltae slightly
% -------------------------------------------------------------------------

ddeltae = deg2rad(1.0);   % small perturbation: +/- 1 deg

% Tail CL at trim +/- elevator perturbation
CLt_plus  = CL0 + at*(alpha_t_trim - it) + CLdeltae*(deltae_trim + ddeltae);
CLt_minus = CL0 + at*(alpha_t_trim - it) + CLdeltae*(deltae_trim - ddeltae);

% Tail induced drag at perturbed elevator values
CDi_t_plus  = (Sth/Sw) * (CLt_plus^2)  / (pi*Art*e);
CDi_t_minus = (Sth/Sw) * (CLt_minus^2) / (pi*Art*e);

% Total drag coefficient at perturbed elevator values
% keep parasite drag and wing induced drag fixed at trim point
CD_plus  = CDp_trim + CDi_w_trim + CDi_t_plus;
CD_minus = CDp_trim + CDi_w_trim + CDi_t_minus;

% central difference derivative
CD_deltae = (CD_plus - CD_minus) / (2*ddeltae);   % [1/rad]

fprintf("CD_deltae (central diff) = %.6f 1/rad\n", CD_deltae);

%% ------------------------------------------------------------------------
% Optional plots for report / sanity check
% -------------------------------------------------------------------------

% Local CD vs alpha fit plot
figure;
plot(rad2deg(alpha(fit_mask)), CD_total(fit_mask), 'bo', 'LineWidth', 1.5); hold on;
plot(rad2deg(alpha(fit_mask)), polyval(p_CD_alpha, alpha(fit_mask)), 'r-', 'LineWidth', 2);
grid on;
xlabel('\alpha [deg]');
ylabel('C_D [-]');
title('Local Linear Fit of C_D vs \alpha at Trim');
legend('Data near trim','Best-fit line','Location','best');

% Plot the elevator perturbation points
figure;
plot(rad2deg([deltae_trim-ddeltae, deltae_trim, deltae_trim+ddeltae]), ...
     [CD_minus, CD_trim, CD_plus], 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('\delta_e [deg]');
ylabel('C_D [-]');
title('Local C_D Change with Elevator Deflection at Trim');
















% ========================================================================
% 12) PLOTS
% ========================================================================
figure;
plot(V, Dp, 'LineWidth', 2); hold on;
plot(V, Di, 'LineWidth', 2);
plot(V, D,  'LineWidth', 2);
plot(V_LDmax, D(idx), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]'); ylabel('Drag [lb]');
title('Parasite vs Induced vs Total Drag vs Velocity');
legend('Dp','Di','D','At (L/D)_{max}','Location','best');

figure;
plot(V, LD, 'LineWidth', 2); hold on;
plot(V_LDmax, LDmax, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]'); ylabel('L/D [-]');
title('L/D vs Velocity');
legend('L/D','(L/D)_{max}','Location','best');

% figure;
% plot(V, Power_required_hp, 'LineWidth', 2); hold on;
% plot(V, Power_available_hp, '--', 'LineWidth', 2);
% grid on;
% xlabel('Velocity V [ft/s]'); ylabel('Power [hp]');
% title('Power Required vs Power Available');
% xlim([30 115])
% ylim([0 3])
% legend('Power Required','Power Available','Location','best');

figure;
plot(V, rad2deg(deltae), 'LineWidth', 2); hold on;
plot(V_LDmax, deltae_LDmax_deg, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]'); ylabel('\delta_e [deg]');
title('Elevator Deflection vs Velocity (Trim from Cm=0)');
legend('\delta_e','At (L/D)_{max}','Location','best');

% figure;
% plot(V, D, 'LineWidth', 2); hold on;
% plot(V, Thrust_available, '--', 'LineWidth', 2);
% grid on;
% xlabel('Velocity V [ft/s]');
% ylabel('Force [lb]');
% title('Drag Required vs Thrust Available');
% legend('Drag Required','Thrust Available','Location','best');


% ========================================================================
% CL^(3/2)/CD vs Velocity (Endurance Metric)
% ========================================================================
figure;
plot(V, CL32_over_CD, 'LineWidth', 2); hold on;
plot(V_endurance, CL32CD_max, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]');
ylabel('C_L^{3/2}/C_D [-]');
title('Endurance Metric vs Velocity');
legend('C_L^{3/2}/C_D','Maximum Endurance Point','Location','best');

% ========================================================================
% CL vs CD (Drag Polar)
% ========================================================================
figure;
plot(CD_total, CL, 'LineWidth', 2);
grid on;
xlabel('C_D [-]');
ylabel('C_L [-]');
title('C_L vs C_D (Drag Polar)');
ylim([0 1.5])


% ========================================================================
% MULTI-RPM POWER AVAILABLE vs POWER REQUIRED
% ========================================================================
figure;
plot(V, Power_required_hp, 'k', 'LineWidth', 3); hold on;

plot(V, Pavail_2000,  'LineWidth', 1.5);
plot(V, Pavail_3000,  'LineWidth', 1.5);
plot(V, Pavail_4000,  'LineWidth', 1.5);
plot(V, Pavail_5000,  'LineWidth', 1.5);
plot(V, Pavail_6000,  'LineWidth', 1.5);
plot(V, Pavail_7000,  'LineWidth', 1.5);
plot(V, Pavail_8000,  'LineWidth', 1.5);
plot(V, Pavail_9000,  'LineWidth', 1.5);
plot(V, Pavail_10000, 'LineWidth', 1.5);



% Vertical markers


if ~isnan(V_top_achievable)
    xline(V_top_achievable, '--r', sprintf('Top Speed = %.1f ft/s', V_top_achievable), 'LineWidth', 1.5);
end

grid on;
xlabel('Velocity V [ft/s]');
ylabel('Power [hp]');
title('Power Required vs Power Available for All RPM Conditions');
legend('Power Required', ...
       '2000 RPM','3000 RPM','4000 RPM','5000 RPM','6000 RPM', ...
       '7000 RPM','8000 RPM','9000 RPM','10000 RPM', ..., ...
       'Location','eastoutside');




% ========================================================================
% PRINT PROP EFFICIENCY AT IMPORTANT SPEEDS
% ========================================================================
V_check = [55 60 78 V_top_achievable];
labels  = {'Stall Speed','Endurance Speed','Max CL/CD Speed','Top Speed'};

eta_2000_chk  = interp1(V, eta_p_2000,  V_check, 'pchip', NaN);
eta_3000_chk  = interp1(V, eta_p_3000,  V_check, 'pchip', NaN);
eta_4000_chk  = interp1(V, eta_p_4000,  V_check, 'pchip', NaN);
eta_5000_chk  = interp1(V, eta_p_5000,  V_check, 'pchip', NaN);
eta_6000_chk  = interp1(V, eta_p_6000,  V_check, 'pchip', NaN);
eta_7000_chk  = interp1(V, eta_p_7000,  V_check, 'pchip', NaN);
eta_8000_chk  = interp1(V, eta_p_8000,  V_check, 'pchip', NaN);
eta_9000_chk  = interp1(V, eta_p_9000,  V_check, 'pchip', NaN);
eta_10000_chk = interp1(V, eta_p_10000, V_check, 'pchip', NaN);

fprintf('\n==================== PROP EFFICIENCY AT KEY SPEEDS ====================\n');
for k = 1:length(V_check)
    fprintf('\n%s at V = %.2f ft/s:\n', labels{k}, V_check(k));
    fprintf('  eta_2000   = %.3f\n', eta_2000_chk(k));
    fprintf('  eta_3000   = %.3f\n', eta_3000_chk(k));
    fprintf('  eta_4000   = %.3f\n', eta_4000_chk(k));
    fprintf('  eta_5000   = %.3f\n', eta_5000_chk(k));
    fprintf('  eta_6000   = %.3f\n', eta_6000_chk(k));
    fprintf('  eta_7000   = %.3f\n', eta_7000_chk(k));
    fprintf('  eta_8000   = %.3f\n', eta_8000_chk(k));
    fprintf('  eta_9000   = %.3f\n', eta_9000_chk(k));
    fprintf('  eta_10000  = %.3f\n', eta_10000_chk(k));
end


% ========================================================================
% 13) TABLE OUTPUT (NO Reynolds; INCLUDE deltae + moments)
% ========================================================================
alpha_deg   = rad2deg(alpha);
alpha_t_deg = rad2deg(alpha_t);
deltae_deg  = rad2deg(deltae);

% Moment breakdown (kept explicit for team)
Cmcgw = Cmacw + aw.*alpha.*(hcg - hacw);                 % wing part about CG
Cmcgt = -(VH*at*(1-downwash)).*alpha + at*VH*it;         % tail part about CG
Cmcg  = Cmcgw + Cmcgt;                                   % total

T = table(V(:), CL(:), alpha_deg(:), alpha_t_deg(:), ...
          deltae_deg(:), CLw(:), CLt(:), ...
          CDp(:), CDi_w(:), CDi_t(:), CD_total(:), ...
          Dp(:), Di(:), D(:), LD(:), ...
          Cmcgw(:), Cmcgt(:), Cmcg(:), ...
          'VariableNames', ...
          {'V_ft_s','CL','alpha_deg','alpha_t_deg', ...
           'deltae_deg','CLw','CLt', ...
           'CDp','CDi_w','CDi_t','CD_total', ...
           'Dp_lb','Di_lb','D_lb','L_over_D', ...
           'Cmcgw','Cmcgt','Cmcg'});

disp(T);

%% ============================ END SCRIPT ================================










%% Stability Derivatves 
%% ========================================================================
%  MAE 154A - STABILITY DERIVATIVES (SIMPLE + MATCHES YOUR VARIABLE NAMES)
%  Paste this near the END of YOUR script (after you computed aw, at, VH_half,
%  lt_half, CL0, CLalpha, Cm0, Cmalpha, CLdeltae, Cmdeltae, etc.)
%
%  Uses YOUR variables:
%   Sw, bw, Cw, e, cla, aw, at, downwash, CL0, CLalpha, Cm0, Cmalpha,
%   CLdeltae, Cmdeltae, VH_half, lt_half, x_cg_half, x_wle, H_fuse, CDp, idx
%
%  VT assumptions requested:
%   - NACA 0012 (handled via thickness in drag; for linear slopes we use cla)
%   - S_vt = 0.08 * Sw
%   - VT has NO taper (lambda_vt = 1.0)
%   - VT trailing edge flush with END of plane at x = 4 ft from nose
%
%  First-pass aero assumptions:
%   eta_t = eta_v = 1, sidewash = 1, no ailerons
% ========================================================================

%% ------------------ 0) SETTINGS (FIRST-PASS) -----------------------------
eta_t    = 1.0;        % HT dynamic pressure ratio
eta_v    = 1.0;        % VT dynamic pressure ratio
sidewash = 1.0;        % (1 + dσ/dβ) ~ 1 first-pass
tau_r    = 0.5;        % rudder effectiveness (like your tau)
x_end    = 4.0;        % [ft] end of plane from nose (YOU SAID 4 ft)

lambda_vt = 1.0;       % no taper
hacv      = 0.25;      % VT AC at 0.25 chord

% Use your design CG case
x_cg = x_cg_half;      % [ft]
VH   = VH_half;        % [-]
lt   = lt_half;        % [ft]
c    = Cw;             % [ft]
b    = bw;             % [ft]
S    = Sw;             % [ft^2]

%% ------------------ 1) VERTICAL TAIL GEOMETRY ----------------------------
% Your requirement: S_vt = 8% of wing area
% (overwrite any older S_vt you had so it matches requirement)
S_vt = 0.08*Sw;

% Use your AR_vt (already in your script)
% (If you forgot to define it, uncomment the next line)
% AR_vt = 1.8;

% VT geometry from S_vt and AR_vt
b_vt  = sqrt(S_vt*AR_vt);                          % [ft] VT span/height
Cr_vt = (2*S_vt)/(b_vt*(1 + lambda_vt));           % [ft] VT root chord

% VT AC location with TE flush at x_end:
% x_vac = x_end - (1 - hacv)*Cr_vt = x_end - 0.75*Cr_vt
x_vac = x_end - (1 - hacv)*Cr_vt;                  % [ft]

% Moment arm CG -> VT AC
lv = x_vac - x_cg;                                 % [ft]

% Vertical tail volume coefficient
Vv = (S_vt*lv)/(Sw*bw);                             % [-]

% VT 3D lift slope (same formula style as your aw/at)
av = cla / (1 + (cla/(pi*e*AR_vt)));               % [1/rad]

% VT vertical offset from CG  sits on top of fuselage)
zv = (H_fuse/2) + (b_vt/2);                        % [ft]

%% ------------------ 2) CD0 FOR Cn_r APPROX (OPTIONAL) --------------------
% Some approximations include -(CD0/4). We'll use parasite drag at your L/Dmax
CD0 = 0;
if exist('CDp','var') && exist('idx','var') && idx>=1 && idx<=length(CDp)
    CD0 = CDp(idx);
end

%% ========================================================================
% 3) LONGITUDINAL DERIVATIVES
% You already compute:
%   CL0, CLalpha, Cm0, Cmalpha, CLdeltae, Cmdeltae
% Here we add: CL_q, CL_adot, Cm_q, Cm_adot
% ========================================================================
CL_q    =  2*eta_t*VH*at;
CL_adot =  2*eta_t*VH*at*downwash;      % uses dε/dα

Cm_q    = -2*eta_t*(lt/c)*VH*at;
Cm_adot = (lt/c)*(-2*eta_t*VH*at*downwash);

%% ========================================================================
% 4) LATERAL-DIRECTIONAL DERIVATIVES (VT-dominated first-pass)
% ========================================================================
CYbeta = -eta_v*(S_vt/Sw)*av*sidewash;
Cnbeta =  eta_v*Vv*av*sidewash;
Clbeta = -eta_v*(zv/bw)*(S_vt/Sw)*av*sidewash;

Clp = -aw/6;

CYr =  2*eta_v*Vv*av;
CYp = -(8/(3*pi))*eta_v*((b_vt*S_vt)/(bw*Sw))*av;

Clr = (CL0/4) + 2*eta_v*(zv/bw)*Vv*av;
Cnr = -(CD0/4) - 2*eta_v*(lv/bw)*Vv*av;

Cnp = 0;  % set 0 first-pass (needs better model)

%% ========================================================================
% 5) CONTROL DERIVATIVES
% No ailerons => set aileron derivatives to 0.
% Rudder => first-pass using VT
% ========================================================================
CYda = 0;  Clda = 0;  Cnda = 0;

CYdr =  eta_v*(S_vt/Sw)*av*tau_r;
Cndr = -eta_v*Vv*av*tau_r;
Cldr = -eta_v*(zv/bw)*(S_vt/Sw)*av*tau_r;

%% ========================================================================
% 6) PRINT CLEAN OUTPUT
% (All variables printed are defined in THIS block or already in your script)
% ========================================================================
fprintf("\n==================== STABILITY DERIVATIVES (MATCHING YOUR CODE) ====================\n");

fprintf("\n--- Longitudinal (you have + new dynamics) ---\n");
fprintf("CL0        = %+ .6f\n", CL0);
fprintf("CLalpha    = %+ .6f\n", CLalpha);
fprintf("CLdeltae   = %+ .6f\n", CLdeltae);
fprintf("Cm0        = %+ .6f\n", Cm0);
fprintf("Cmalpha    = %+ .6f\n", Cmalpha);
fprintf("Cmdeltae   = %+ .6f\n", Cmdeltae);

fprintf("CL_q       = %+ .6f\n", CL_q);
fprintf("CL_adot    = %+ .6f\n", CL_adot);
fprintf("Cm_q       = %+ .6f\n", Cm_q);
fprintf("Cm_adot    = %+ .6f\n", Cm_adot);

fprintf("\n--- Lateral/Directional (first-pass VT dominated) ---\n");
fprintf("CYbeta     = %+ .6f\n", CYbeta);
fprintf("Clbeta     = %+ .6f\n", Clbeta);
fprintf("Cnbeta     = %+ .6f\n", Cnbeta);

fprintf("Clp        = %+ .6f\n", Clp);
fprintf("CYp        = %+ .6f\n", CYp);
fprintf("CYr        = %+ .6f\n", CYr);
fprintf("Clr        = %+ .6f\n", Clr);
fprintf("Cnr        = %+ .6f   (used CD0 = %.6f)\n", Cnr, CD0);
fprintf("Cnp        = %+ .6f   (set 0 first-pass)\n", Cnp);

fprintf("\n--- Controls (no ailerons; rudder included) ---\n");
fprintf("CYda       = %+ .6f (no ailerons)\n", CYda);
fprintf("Clda       = %+ .6f (no ailerons)\n", Clda);
fprintf("Cnda       = %+ .6f (no ailerons)\n", Cnda);

fprintf("CYdr       = %+ .6f\n", CYdr);
fprintf("Cldr       = %+ .6f\n", Cldr);
fprintf("Cndr       = %+ .6f\n", Cndr);

%% ========================================================================
% 7) WHAT YOU ARE STILL MISSING (WITH YOUR CURRENT INFO)
% ========================================================================
fprintf("\n==================== STILL MISSING / NOT WELL-DEFINED YET ====================\n");
fprintf("1) Aileron derivatives: CYda, Clda, Cnda (we set them to 0 since you said no ailerons).\n");
fprintf("2) Better Cnp (we set to 0 first-pass; needs more detailed lateral model).\n");
fprintf("3) Xu, Zu, Mu (speed derivatives) unless you model thrust vs speed and CD vs speed at trim.\n");
fprintf("4) Inertias Ix, Iy, Iz, Ixz (needed for full state-space / simulator dynamics).\n");
fprintf("5) Exact sign convention check vs your lecture table (depends on whether your class uses CZ or CL sign).\n");