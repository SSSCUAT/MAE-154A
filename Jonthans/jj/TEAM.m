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
W_total     = 19;             % [lb]
Power = 2.8;            % [hp] power available (flat line first-pass)

% --- Wing geometry ---
Arw     = 7.4;          % [-]
bw      = 5.617;        % [ft]
Sw      = 4.264931025;  % [ft^2]
lambdaw = 0.40;         % [-]

% --- Horizontal tail geometry ---
Sth      = Sw*0.20;     % [ft^2]
Art      = 4.0;         % [-]
bth      = sqrt(Sth*Art); % [ft]
lambdath = 0.50;        % [-]

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
e        = 1;        % [-]
cla      = 5.73;        % [1/rad]
downwash = 0.25;        % [-] dε/dα
tau      = 0.5;         % [-] elevator effectiveness
Cmacw    = -0.05;       % [-] wing Cm about wing AC 

% --- FIXED TAIL INCIDENCE (FIRST PASS) ---
it = 3*pi/180;          % [rad] 

% --- Locations (from nose) ---
x_wle     = 1.0;        % [ft] wing leading edge location
x_cg_dry  = 1.5339;     % [ft]
x_cg_half = 1.5103;     % [ft] 
x_cg_0    = 1.4496;     % [ft]

% --- Fuselage (box: 8in x 8in x 5ft) ---
L_fuse = 5.0;           % [ft]
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
fprintf("SM_dry  = %.4f (%.1f%% MAC)\n", SM_dry,  100*SM_dry);
fprintf("SM_half = %.4f (%.1f%% MAC)\n", SM_half, 100*SM_half);
fprintf("SM_0    = %.4f (%.1f%% MAC)\n", SM_0,    100*SM_0);

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

% Tail lift coefficient using YOUR professor-style form:
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
Power_required_hp = (D.*V) / 550;
Power_available_hp = Power * ones(size(V));

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

figure;
plot(V, Power_required_hp, 'LineWidth', 2); hold on;
plot(V, Power_available_hp, '--', 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]'); ylabel('Power [hp]');
title('Power Required vs Power Available');
legend('Power Required','Power Available','Location','best');

figure;
plot(V, rad2deg(deltae), 'LineWidth', 2); hold on;
plot(V_LDmax, deltae_LDmax_deg, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]'); ylabel('\delta_e [deg]');
title('Elevator Deflection vs Velocity (Trim from Cm=0)');
legend('\delta_e','At (L/D)_{max}','Location','best');

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

% VT vertical offset from CG (first-pass: sits on top of fuselage)
% You already have H_fuse in your script.
zv = (H_fuse/2) + (b_vt/2);                        % [ft]

%% ------------------ 2) CD0 FOR Cn_r APPROX (OPTIONAL) --------------------
% Some approximations include -(CD0/4). We'll use parasite drag at your L/Dmax
% point (idx) as "CD0" (first-pass).
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