%% 154A Project
clear; clc;

% Aerodynamic Variables
CL = 1.5;        % CL
AR = 7.4;
CD_p = 0.027;
e = 0.80;
S = 9.3;         % Cessna 182 wing area (ft^2) Miniture
W = 50;          %lbs 
Power = 2.8  ;   %Power Motor in HP
p = Power*550 ;  %Power in lb(ft/s^2)

% Induced Drag
CDi = (CL^2) / (pi * e * AR);

%% Inputs
v = 55;         % ft/s 
L = 0.767 * 3.28084;  % chord in ft, found online

%% Sea Level Constants
rho = 0.0023769;    % slug/ft^3
mu  = 3.737e-7;     % slug/(ft·s)

%% Reynolds Number
Re = (rho * v * L) / mu;

%% Power Required
q = 0.5 * rho * v^2;
cd_t = CDi + CD_p;
D = q * S * cd_t;        % lbf
P = D * v;               % ft*lbf/s
P_hp = P / 550;          % horsepower

%% Output
fprintf('Velocity: %.1f ft/s\n', v);
fprintf('Reynolds Number: %.2e\n', Re);
fprintf('Power Required: %.1f hp\n', P_hp);
fprintf('CDi: %.4f\n', CDi);

%% Lift Drag Curve 
% Drag vs Velocity Curve
v1 = linspace(30,200,400);   % ft/s

CL1  = (2*W) ./ (rho .* v1.^2 .* S);
CDi1 = (CL1.^2) ./ (pi * e * AR);

q1 = 0.5 .* rho .* v1.^2;

D_p = q1 .* S .* CD_p;     % parasite drag
D_I = q1 .* S .* CDi1;    % induced drag
D1  = D_p + D_I;          % total drag
Thurst = p./(v1);

figure 
hold on
plot(v1, D1,  'b-',  'LineWidth', 2, 'DisplayName', 'Total Drag')
plot(v1, D_p, 'r--', 'LineWidth', 2, 'DisplayName', 'Parasite Drag')
plot(v1, D_I, 'g:',  'LineWidth', 2, 'DisplayName', 'Induced Drag')
%plot(v1, Thurst, 'white', 'LineWidth',2, 'DisplayName','Thurst')
hold off

xlabel('Velocity (ft/s)')
ylabel('Drag (lbf)')
title('Drag vs Velocity')
legend('Location','best')
grid on

%% power curve 
v_1 = linspace(5,200,400);
CL2  = (2*W) ./ (rho .* v_1.^2 .* S);
CDi2 = (CL2.^2) ./ (pi * e * AR);

q2 = 0.5 .* rho .* v_1.^2;

D_p2 = q2 .* S .* CD_p;     % parasite drag
D_I2 = q2 .* S .* CDi2;    % induced drag
D2  = D_p2 + D_I2;          % total drag



Power_Required1 = D2 .* v_1;          % ft*lbf/s
Power_Required_hp = Power_Required1 / 550;

P_available_hp = Power * ones(size(v_1));  % turn Power avalaible in a vector 

figure
hold on
plot(v_1, Power_Required_hp, 'b-', 'LineWidth', 2, ...
     'DisplayName', 'Power Required')
plot(v_1, P_available_hp, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'Power Available')
hold off

xlabel('Velocity (ft/s)')
ylabel('Power (hp)')
title('Power Required vs Velocity')
legend('Location','best')
grid on
