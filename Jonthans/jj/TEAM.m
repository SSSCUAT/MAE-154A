%% ========================================================================
clear 
clc 
clear all 
% 0) ATMOSPHERE (SEA LEVEL ASSUMPTION)
% ========================================================================
rho = 0.0023769;        % Air density at sea level [slug/ft^3]
mu  = 3.737e-7;         % Dynamic viscosity at sea level [slug/(ft*s)]
a   = 1116.45;          % Speed of sound at sea level [ft/s]

%% ========================================================================
% 1) AIRCRAFT INPUTS (EDIT THESE WHEN DESIGN CHANGES)
% ========================================================================

% -----------------------------
% Weight (used for drag sweep)
% -----------------------------
W   = 25;                     % Aircraft weight [lb] (update as weight estimate changes)

% -----------------------------
% Wing geometry (KNOWN)
% -----------------------------
Arw     = 7.4;                % Wing aspect ratio [-]
bw      = 5.617;              % Wing span [ft]
Sw      = 4.264931025;        % Wing planform area [ft^2]
lambdaw = 0.40;               % Wing taper ratio ct/cr [-]

% -----------------------------
% Horizontal tail geometry (KNOWN)
% -----------------------------
Sth      = 0.5331163781;      % Horizontal tail area [ft^2]
bth      = 1.460296378;       % Horizontal tail span [ft]
lambdath = 0.50;              % Tail taper ratio ct/cr [-]
Art      = 4.0;               % Tail aspect ratio [-] (given/assumed)

% -----------------------------
% Aero assumptions (EDIT AS NEEDED)
% -----------------------------
e        = 0.85;              % Oswald efficiency factor [-] (guess)
cla      = 5.73;              % 2D airfoil lift curve slope [1/rad] (CHECK units)
downwash = 0.25;              % Downwash gradient dε/dα [-]
it       = -1*pi/180;         % Tail incidence angle [rad] (initial guess)
Cmacw    = -0.05;             % Wing Cm about wing AC (approx from airfoil data)
tau      = 0.5;               % Elevator effectiveness τe [-] (guess) (used only if you later trim with δe)
CL_MAX = 1.5
% -----------------------------
% Layout / CG reference positions (from CG build-up)
% -----------------------------
x_wle    = 0.75;              % Wing LE station from nose [ft]
x_cg_str = 1.7147;            % Structural CG station from nose [ft]
x_cg_tot = 1.4103;            % CG with avionics station from nose [ft]
x_cg_0   = 1.5071;            % CG with payload+fuel station from nose [ft]

% -----------------------------
% Fuselage geometry (for drag)
% Given: 8in x 8in x 5ft box
% -----------------------------
L_fuse = 5.0;                 % Fuselage length [ft]
W_fuse = 8/12;                % Fuselage width [ft]
H_fuse = 8/12;                % Fuselage height [ft]
Power = 2.5;























% Estimated fuselage wetted area [ft^2] (rectangular prism):
% Swet_f = 2(LW + LH + WH)
Swet_f = 2*(L_fuse*W_fuse + L_fuse*H_fuse + W_fuse*H_fuse);

% Equivalent "diameter" used in fuselage fineness ratio for Raymer form factor.
% For a square-ish fuselage, using width is a reasonable approximation.
d_fuse = W_fuse;              % [ft]

%% ========================================================================
% 2) PLANFORM GEOMETRY (ROOT, TIP, MAC) - WING + TAIL
% ========================================================================

% ---- Wing chords (trapezoid geometry) ----
Crw = (2*Sw)/(bw*(1 + lambdaw));                 % Wing root chord [ft]
Ctw = lambdaw * Crw;                             % Wing tip chord [ft]
Cw  = (2/3) * Crw * ((1 + lambdaw + lambdaw^2)/(1 + lambdaw));  % Wing MAC [ft]

% ---- Tail chords (trapezoid geometry) ----
Crth = (2*Sth)/(bth*(1 + lambdath));             % Tail root chord [ft]
Ctth = lambdath * Crth;                          % Tail tip chord [ft]
Cth  = (2/3) * Crth * ((1 + lambdath + lambdath^2)/(1 + lambdath)); % Tail MAC [ft]

fprintf("\n==================== GEOMETRY ====================\n");
fprintf("Wing:  Cr=%.4f ft, Ct=%.4f ft, MAC Cw=%.4f ft\n", Crw, Ctw, Cw);
fprintf("Tail:  Cr=%.4f ft, Ct=%.4f ft, MAC Cth=%.4f ft\n", Crth, Ctth, Cth);
fprintf("Fuse:  L=%.2f ft, W=%.3f ft, H=%.3f ft, Swet_f=%.3f ft^2\n", L_fuse, W_fuse, H_fuse, Swet_f);

%% ========================================================================
% 3) 3D LIFT CURVE SLOPES (FINITE WING CORRECTION)
% ========================================================================
% 3D lift slope formula:
% a = a0 / (1 + a0/(pi*e*AR))
aw = cla / (1 + (cla/(pi*e*Arw)));   % Wing 3D lift curve slope [1/rad]
at = cla / (1 + (cla/(pi*e*Art)));   % Tail 3D lift curve slope [1/rad]
%Just fyi effeicny factor is a guess 
fprintf("\n==================== LIFT SLOPES ====================\n");
fprintf("aw (wing 3D) = %.4f 1/rad\n", aw);
fprintf("at (tail 3D) = %.4f 1/rad\n", at);

%% ========================================================================
% 4) LONGITUDINAL REFERENCE LOCATIONS (MEASURED FROM WING LE)
% ========================================================================
% Wing aerodynamic center location (fraction of wing MAC)
hacw = 0.25;                      % Wing AC at 25% MAC (standard subsonic assumption)

% Wing AC from wing LE (dimensional)
x_ac_w_LE = hacw*Cw;              % [ft]

% Tail arm assumption:
% Lh is distance from wing AC to tail AC (design/layout choice)
% Here: assume Lh = 4*Cw (update once layout is finalized)
Lh = 4*Cw;                        % [ft] wing AC -> tail AC

% Tail AC from wing LE (dimensional)
x_ac_t_LE = x_ac_w_LE + Lh;       % [ft]

% Tail AC nondimensional from wing LE
h_ac_t = x_ac_t_LE / Cw;          % [-]

%% ========================================================================
% 5) NEUTRAL POINT (hn) + STATIC MARGIN (SM)
% ========================================================================
% Neutral point for wing + tail:
hn = (hacw + h_ac_t*((Sth/Sw)*(at/aw))*(1-downwash)) / (1 + (Sth/Sw)*(at/aw)*(1-downwash));

% CG locations nondimensional from wing LE:
hcg_str = (x_cg_str - x_wle)/Cw;
hcg_tot = (x_cg_tot - x_wle)/Cw;
hcg_0   = (x_cg_0   - x_wle)/Cw;

% Static margin: SM = hn - hcg (positive = stable)
SM_str = hn - hcg_str;
SM_tot = hn - hcg_tot;
SM_0   = hn - hcg_0;

fprintf("\n==================== STABILITY ====================\n");
fprintf("Neutral point hn = %.4f (fraction of MAC from wing LE)\n", hn);
fprintf("hcg_str = %.4f  --> SM_str = %.4f (%.1f%% MAC)\n", hcg_str, SM_str, 100*SM_str);
fprintf("hcg_tot = %.4f  --> SM_tot = %.4f (%.1f%% MAC)\n", hcg_tot, SM_tot, 100*SM_tot);
fprintf("hcg_0   = %.4f  --> SM_0   = %.4f (%.1f%% MAC)\n", hcg_0,   SM_0,   100*SM_0);

%% ========================================================================
% 6) TAIL ARM lt (CG -> TAIL AC) AND TAIL VOLUME VH
% ========================================================================
% Convert wing/tail AC to nose stations, then compute tail arm lt = x_tac - x_cg.
x_wac = x_wle + x_ac_w_LE;        % Wing AC station from nose [ft]
x_tac = x_wac + Lh;               % Tail AC station from nose [ft]

lt_str = x_tac - x_cg_str;        % Tail arm for structural CG [ft]
lt_tot = x_tac - x_cg_tot;        % Tail arm for avionics CG [ft]
lt_0   = x_tac - x_cg_0;          % Tail arm for payload+fuel CG [ft]

% Tail volume coefficient: VH = (St * lt) / (Sw * Cw)
VH_str = (Sth*lt_str)/(Sw*Cw);
VH_tot = (Sth*lt_tot)/(Sw*Cw);
VH_0   = (Sth*lt_0)  /(Sw*Cw);

fprintf("\n==================== TAIL VOLUME ====================\n");
fprintf("Tail AC station x_tac = %.4f ft\n", x_tac);
fprintf("lt_str = %.4f ft --> VH_str = %.4f\n", lt_str, VH_str);
fprintf("lt_tot = %.4f ft --> VH_tot = %.4f\n", lt_tot, VH_tot);
fprintf("lt_0   = %.4f ft --> VH_0   = %.4f\n", lt_0,   VH_0);

%% ========================================================================
% 7) LONGITUDINAL COEFFICIENTS (ALL OF THEM)
% ========================================================================

% Choose loading case for trim/stability calculations
VH  = VH_tot;        % Tail volume coefficient
hcg = hcg_tot;       % CG location (nondim from wing LE)

% Given/assumed
% downwash = dε/dα (we use 0.25 from Cessna-type data)
% Cmacw = wing Cm about wing AC (airfoil data)
% tau = elevator effectiveness (guess)

% --- Lift coefficients ---
CLalpha   = aw + at*(Sth/Sw)*(1 - downwash);   % dCL/dalpha (aircraft)
CL0       = -at*(Sth/Sw)*it;                   % CL at alpha = 0 (tail incidence)

% --- Moment coefficients (about CG) ---
Cmalpha   = -CLalpha*(hn - hcg);               % dCm/dalpha (stability)
Cm0       = Cmacw + VH*at*it;                  % Cm at alpha = 0

% --- Moment breakdown pieces (functions of alpha; compute later in sweep) ---
 %Cmcgw = Cmacw + aw*alpha*(hcg - hacw);        % wing moment about CG
 %Cmcgt = -(VH*at*(1-downwash))*alpha + at*VH*it; % tail moment about CG
 %Cmcg  = Cmcgw + Cmcgt;                        % total moment about CG

% --- Elevator trim equation (deltae is a function of CL) ---
% deltae = -((Cm0*CLalpha) + (Cmalpha*CL)) ./ ((CLalpha*Cmdeltae) - (Cmalpha*CLdeltae));
% Control derivatives (if/when you use elevator):
CLdeltae = at*(Sth/Sw)*tau;     % dCL/dδe (usually +)
Cmdeltae = -at*VH*tau;          % dCm/dδe (usually -)

% NOTE: To actually compute δe at a flight condition, you must define CL and use:
% deltae = -((Cm0*CLalpha) + (Cmalpha*CL)) / ((CLalpha*Cmdeltae) - (Cmalpha*CLdeltae));

%% ========================================================================
% 8) VELOCITY SWEEP (DRAG + L/D)  [THIS IS THE MAIN PERFORMANCE PART]
% ========================================================================
% We sweep velocity, compute required CL from L=W at each speed, then compute:
%   - Trim alpha and tail alpha_t from your linear model
%   - CDp from Raymer build-up (Re-dependent)
%   - CDi from induced drag using CLw and CLt
%   - Drag forces and L/D

% -----------------------------
% 8.1 Velocity range
% -----------------------------
V = linspace(30, 300, 500);      % [ft/s] (edit as desired)
q = 0.5*rho.*V.^2;               % Dynamic pressure [lb/ft^2]
M = V./a;                        % Mach number [-]

% -----------------------------
% 8.2 Required CL for level flight
% -----------------------------
CL = W ./ (q.*Sw);               % Total aircraft CL required for L=W

% -----------------------------
% 8.3 Solve for alpha, tail AoA, and CLw/CLt (same model you've been using)
% -----------------------------
alpha = CL ./ aw;          % wing-only alpha (rad)
alpha_t = (1 - downwash).*alpha + it;

CLw = aw .* alpha;
CLt = at .* alpha_t;

% -----------------------------
% 8.4 Reynolds numbers (for skin friction)
% -----------------------------
Re_w = rho.*V.*Cw   ./ mu;       % Re based on wing MAC
Re_t = rho.*V.*Cth  ./ mu;       % Re based on tail MAC
Re_f = rho.*V.*L_fuse ./ mu;     % Re based on fuselage length

% -----------------------------
% 8.5 Parasite drag build-up (Raymer style)
% -----------------------------
% Wetted areas (simple first-pass approximations)
Swet_w = 2*Sw;                   % wing wetted area ~ 2*planform
Swet_t = 2*Sth;                  % tail wetted area ~ 2*planform
% Swet_f already computed from box geometry

% Skin friction coefficients (Raymer Eq. style)
Cf_w = 0.455 ./ ((log10(Re_w)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_t = 0.455 ./ ((log10(Re_t)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_f = 0.455 ./ ((log10(Re_f)).^2.58 .* (1 + 0.144*M.^2).^0.65);

% Form factor settings (surface thickness effects, sweep effects)
tc  = 0.12;                      % thickness-to-chord ratio t/c [-] (guess)
xcm = 0.30;                      % location of max thickness x/c [-] (guess)
Lam = 0;                         % sweep at max thickness [rad] (0 for unswept)
K_surf = (1 + (0.6/xcm)*tc + 100*tc^4) .* (1.34*M.^0.18 .* (cos(Lam)).^0.28);

% Interference factors (how junctions increase drag)
Q_w = 1.0;                       % wing interference factor
Q_t = 1.05;                      % tail interference factor
Q_f = 1.0;                       % fuselage factor

% Fuselage form factor (Raymer smooth-body approximation)
f = L_fuse/d_fuse;               % fineness ratio
K_f = 0.9 + 5/f^1.5 + f/400;

% Misc drag (gaps, protuberances) and leakage/protuberance allowances
CD_misc = 0.002;
CD_LP   = 0.002;

% Total parasite drag coefficient CDp
CDp = ...
    K_surf.*Q_w.*Cf_w.*(Swet_w/Sw) + ...
    K_surf.*Q_t.*Cf_t.*(Swet_t/Sw) + ...
    K_f.*Q_f.*Cf_f.*(Swet_f/Sw) + ...
    CD_misc + CD_LP;

% -----------------------------
% 8.6 Induced drag (wing + tail)
% -----------------------------
CDi_w = (CLw.^2) ./ (pi*Arw*e);               % wing induced drag
CDi_t = (Sth/Sw) .* (CLt.^2) ./ (pi*Art*e);   % tail induced drag (scaled by St/Sw)

CDi_total = CDi_w + CDi_t;
CD_total  = CDp + CDi_total;

% -----------------------------
% 8.7 Forces and performance metrics
% -----------------------------
Dp = q .* Sw .* CDp;            % Parasite drag force [lb]
Di = q .* Sw .* CDi_total;      % Induced drag force [lb]
D  = q .* Sw .* CD_total;       % Total drag force [lb]

L  = q .* Sw .* CL;             % Lift force [lb] (should equal W)
LD = L ./ D;                    % Lift-to-drag ratio [-]
%% Power curve
Power_Required1 = D .* V;          % ft*lbf/s
Power_Required_hp = Power_Required1 / 550;
mask_stall = CL > CL;
%Power_Required_hp(mask_stall) = NaN;% Dont plot points for velocity
%understall velocity(past max cl value

P_available_hp = Power * ones(size(V));  % turn Power avalaible in a vector 

figure
hold on
plot(V, Power_Required_hp, 'b-', 'LineWidth', 2, ...
     'DisplayName', 'Power Required')
plot(V, P_available_hp, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'Power Available')
hold off

xlabel('Velocity (ft/s)')
ylabel('Power (hp)')
title('Power Required vs Velocity')
legend('Location','best')
grid on
 
%%
% Find (L/D)_max and corresponding velocity
[LDmax, idx] = max(LD);
V_LDmax = V(idx);

fprintf('\n==================== VELOCITY SWEEP RESULTS ====================\n');
fprintf('(L/D)max = %.3f\n', LDmax);
fprintf('V at (L/D)max = %.2f ft/s\n', V_LDmax);
fprintf('CL at (L/D)max = %.3f\n', CL(idx));
fprintf('CDp at (L/D)max = %.4f,  CDi_total at (L/D)max = %.4f\n', CDp(idx), CDi_total(idx));

%% ========================================================================
% 9) PLOTS (DRAG COMPONENTS AND L/D)
% ========================================================================

figure;
plot(V, Dp, 'LineWidth', 2); hold on;
plot(V, Di, 'LineWidth', 2);
plot(V, D,  'LineWidth', 2);
plot(V_LDmax, D(idx), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]');
ylabel('Drag Force [lb]');
title('Parasite vs Induced vs Total Drag vs Velocity');
legend('Parasite Drag D_p','Induced Drag D_i','Total Drag D','At (L/D)_{max}','Location','best');

figure;
plot(V, LD, 'LineWidth', 2); hold on;
plot(V_LDmax, LDmax, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Velocity V [ft/s]');
ylabel('L/D [-]');
title('L/D vs Velocity');
legend('L/D','(L/D)_{max}','Location','best');

%% ========================================================================
% 10) TABLE OUTPUT (TEAM-FRIENDLY)
% ========================================================================
% This table is useful for copying into reports or checking values at a speed.
alpha_deg   = rad2deg(alpha);
alpha_t_deg = rad2deg(alpha_t);
T = table(V(:), CL(:), alpha_deg(:), alpha_t_deg(:), CLw(:), CLt(:), ...
          CDp(:), CDi_w(:), CDi_t(:), CD_total(:), ...
          Dp(:), Di(:), D(:), LD(:), ...
          Re_w(:), Re_t(:), Re_f(:), ...
          'VariableNames', ...
          {'V_ft_s','CL','alpha_deg','alpha_t_deg','CLw','CLt', ...
           'CDp','CDi_w','CDi_t','CD_total', ...
           'Dp_lb','Di_lb','D_lb','L_over_D', ...
           'Re_w','Re_t','Re_f'});

disp(T);

% ============================ END SCRIPT =================================