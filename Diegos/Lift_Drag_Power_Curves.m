<<<<<<< HEAD
%% 154A Project
clear; clc;

% Aerodynamic Variables
CL = 1.5;        % CL max 
AR = 7.4;
CD_p = 0.027;
e = 0.80;
S = 5.6;         % Cessna 182 wing area (ft^2) Miniture
W = 30;          %lbs 
Power = 2.8  ;   %Power Motor in HP
p = Power*550 ;  %Power in lb(ft/s^2)

% Induced Drag
CDi = (CL^2) / (pi * e * AR);
display(CDi)
%% Inputs
v = 55;         % ft/s 
L = 0.767 * 3.28084;  % chord in ft, found online

%% Sea Level Constants
rho = 0.0023769;    % slug/ft^3
mu  = 3.737e-7;     % slug/(ft·s)

%% Reynolds Number
Re = (rho * v * L) / mu;
=======
>>>>>>> f471e6bbb995297936b74c5b570055cc1b7c9fda

%% Power Required
q = 0.5 * rho * v^2;
cd_t = CDi + CD_p;
D = q * S * cd_t;        % lbf
P = D * v;               % ft*lbf/s
P_hp = P / 550;          % horsepower




%% power curve 
v_1 = linspace(4,205,400);
CL2  = (2*W) ./ (rho .* v_1.^2 .* S);
CDi2 = (CL2.^2) ./ (pi * e * AR);

q2 = 0.5 .* rho .* v_1.^2;

D_p2 = q2 .* S .* CD_p;     % parasite drag
D_I2 = q2 .* S .* CDi2;    % induced drag
D2  = D_p2 + D_I2;          % total drag



Power_Required1 = D2 .* v_1;          % ft*lbf/s
Power_Required_hp = Power_Required1 / 550;
%
mask_stall = CL2 > CL;
%Power_Required_hp(mask_stall) = NaN;% Dont plot points for velocity
%understall velocity(past max cl value

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
 
%% Fuel consumption 
Flight_Time_Min = 45;  % How many minutes you want to fly

% Consumption rate for a 20cc engine at 2.8HP is approx 0.6 oz/min
Fuel_Rate = 0.60

% Calculation
Total_Gas_Needed = Flight_Time_Min * Fuel_Rate;
Safety_Reserve = Total_Gas_Needed * 0.20; % 20% reserve is standard
Total_Tank_Size = Total_Gas_Needed + Safety_Reserve;

% --- Results ---
fprintf('Flight Duration: %d minutes\n', Flight_Time_Min);
fprintf('Gasoline needed: %.1f fl oz\n', Total_Gas_Needed);
fprintf('Recommended tank size (with reserve): %.1f fl oz\n', Total_Tank_Size);