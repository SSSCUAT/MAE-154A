%% MAE 154A - Drag Table vs CL (1.5 down to 0.5)
% -------------------------------------------------------------------------
% We sweep CL from 1.5 down to 0.5
% Then compute velocity from L = W
% Everything else follows from CL
% -------------------------------------------------------------------------

clear; clc;

%% ========================================================================
% 1) SEA LEVEL CONDITIONS
% ========================================================================

rho = 0.0023769;        % Sea level density [slug/ft^3]
mu  = 3.737e-7;         % Sea level viscosity [slug/(ft*s)]
a   = 1116.45;          % Speed of sound [ft/s]

%% ========================================================================
% 2) AIRCRAFT PARAMETERS
% ========================================================================

W   = 25;               % Weight [lb]
Sw  = 4.264931025;      % Wing area [ft^2]
Arw = 7.4;              % Wing aspect ratio [-]
Sth = 0.5331163781;     % Tail area [ft^2]
Art = 4.0;              % Tail aspect ratio [-]

Cw  = 0.80;             % Wing MAC for Reynolds [ft]
Cth = 0.35;             % Tail chord for Reynolds [ft]

e        = 0.85;        % Oswald efficiency [-]
downwash = 0.25;        % dε/dα [-]
it       = -1*pi/180;   % Tail incidence [rad]

cla = 5.73;                             % 2D lift slope [1/rad]
aw  = cla / (1 + (cla/(pi*e*Arw)));      % Wing 3D lift slope [1/rad]
at  = cla / (1 + (cla/(pi*e*Art)));      % Tail 3D lift slope [1/rad]

%% ========================================================================
% 3) CL SWEEP (FROM 1.5 DOWN TO 0.5)
% ========================================================================

N  = 25;                            % Number of rows
CL = linspace(1.5, 0.5, N);         % CL vector

%% ========================================================================
% 4) VELOCITY FROM L = W
% ========================================================================

V = sqrt((2*W) ./ (rho*Sw.*CL));     % Velocity [ft/s]
q = 0.5*rho.*V.^2;                   % Dynamic pressure [lb/ft^2]
M = V./a;                            % Mach number [-]

%% ========================================================================
% 5) SOLVE FOR alpha AND alpha_t
% ========================================================================

alpha = ( CL - (Sth/Sw)*at*it ) ./ ( aw + (Sth/Sw)*at*(1 - downwash) );

alpha_t = (1 - downwash).*alpha + it;

CLw = aw .* alpha;
CLt = at .* alpha_t;

%% ========================================================================
% 6) REYNOLDS NUMBERS
% ========================================================================

Re_w = rho.*V.*Cw ./ mu;
Re_t = rho.*V.*Cth ./ mu;
Re_f = rho.*V.*4.0 ./ mu;   % Using 4 ft fuselage length

%% ========================================================================
% 7) PARASITE DRAG (Raymer Build-Up)
% ========================================================================

Swet_w = 2*Sw;
Swet_t = 2*Sth;
Swet_f = 6.0;

tc = 0.12;
xcm = 0.30;

Lam = 0;
Q_w = 1.0;
Q_t = 1.05;
Q_f = 1.0;

Cf_w = 0.455 ./ ((log10(Re_w)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_t = 0.455 ./ ((log10(Re_t)).^2.58 .* (1 + 0.144*M.^2).^0.65);
Cf_f = 0.455 ./ ((log10(Re_f)).^2.58 .* (1 + 0.144*M.^2).^0.65);

K = (1 + (0.6/xcm)*tc + 100*tc^4) .* (1.34*M.^0.18 .* (cos(Lam)).^0.28);

f = 4.0/0.4;                         % fuselage fineness ratio
K_f = 0.9 + 5/f^1.5 + f/400;

CDp = ...
    K.*Q_w.*Cf_w.*(Swet_w/Sw) + ...
    K.*Q_t.*Cf_t.*(Swet_t/Sw) + ...
    K_f.*Q_f.*Cf_f.*(Swet_f/Sw) + ...
    0.002 + 0.002;

%% ========================================================================
% 8) INDUCED DRAG
% ========================================================================

CDi_w = (CLw.^2) ./ (pi*Arw*e);
CDi_t = (Sth/Sw) .* (CLt.^2) ./ (pi*Art*e);

CD_total = CDp + CDi_w + CDi_t;

D = q .* Sw .* CD_total;

%% ========================================================================
% 9) TABLE OUTPUT
% ========================================================================

T = table(V(:), CL(:), alpha(:), CLw(:), alpha_t(:), ...
          CDp(:), CDi_w(:), CDi_t(:), CD_total(:), D(:), ...
          Re_w(:), Re_t(:), Re_f(:), ...
          'VariableNames', ...
          {'V_ft_s','CL','alpha_rad','CLw','alpha_t_rad',...
           'CDp','CDi_w','CDi_t','CD_total','D_lb',...
           'Re_w','Re_t','Re_f'});

disp(T);