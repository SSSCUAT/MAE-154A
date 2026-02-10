% All inputs for the Pioneer UAV dynamics should be done here. AKA all
% constant block values need to be declared here.

whaleSeed_E = sum(clock)+2;
whaleSeed_N = sum(clock)+3;


weight = 420;     % [lbsf]        Bray pg 33
I_xx = 34.832;    % [slug-ft^2] IAI, MSS notes
I_yy = 67.08;     % [slug-ft^2] IAI, MSS notes
I_zz = 82.22;     % [slug-ft^2] IAI, MSS notes
I_xz = -4.902;    % [slug-ft^2] IAI, MSS notes

bw   =  16.90;    % [ft]        Bray pg 31
cbar =  1.80;     % [ft]        Bray pg 31
Sw   =  30.42;    % [ft^2]      Bray pg 31

thrust_bx = 89.8;     % [lb] - body x thrust  from Bray pg 78/80 thrust avail for 65.3 knots at sea lvl
He=0;   % Engine moment

% gravity acceleration
grav = 32.17405;

mass = weight/grav; % [lbm]



% density
rho = 0.002379;     

% need to initialize 6DOF block and aero coef data
%%%
%initial position in inertial axes [Xe, Ye, Ze]
Posini = [0 0 -1000];

%initial velocity in body axes [U, v, w] - relative to wind i think?
Velini = [118.1 0 0];    % 65.3 knots  - zero elevator trim from Bray pg 35

%initial Euler orientation [roll, pitch, yaw]
% Eulini = [0 6.2*pi/180 0];

Eulini = [0 0 0];

%initial angular velocities [p,q,r]
Angvelini = [0 0 0];

%initial aero angles
alpha0 = 0; %6.2*pi/180;     % 6.2 degrees trim alpha for de = 0  from Bray pg 35
beta0 = 0;

%mass is above
%inertia matrix
Inertia = [I_xx 0 I_xz; 0 I_yy 0; I_xz 0 I_zz];

%%%


