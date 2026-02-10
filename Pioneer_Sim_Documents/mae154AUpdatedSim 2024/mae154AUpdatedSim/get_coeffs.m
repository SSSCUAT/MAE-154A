function out = get_coeffs(in)

alpha = in(1);
alpha_dot = in(2);
beta = in(3);
elevator = in(4);
aileron = in(5);
p = in(6);
q = in(7);
r = in(8);
rudder = in(9);
Vt = in(10);

cbar =  1.80;     % [ft]        Bray pg 31
bw   =  16.90;    % [ft]        Bray pg 31


% aero data - linear stability coefficients
CL0     =    0.385;              % []          Bray pg 33
CL_a    =    4.78;               % [/rad]      Bray pg 33
CL_adot =    2.42;               % [/rad]      Bray pg 33
CL_q    =    8.05;               % [/rad]      Bray pg 33
CL_de   =    0.401;              % [/rad]      Bray pg 33
%CLfa CLfa.dat 0 1               %             Bray pg 50, Table 4.7
%CLfade CLfade.dat 0 1 1         %             Bray pg 41, Table 4.4

CD0     =    0.060;              % []          Bray pg 33
CD_a    =    0.430;              % [/rad]      Bray pg 33
CD_de   =    0.0180;             % [/rad]      Bray pg 33
%CDfa Cdfa.dat 0 1               %             Bray pg 50, Table 4.7
%CDfade CDfade.dat 0 1 1         %             Bray pg 45, Table 4.6

Cm0     =   0.194;              % []          Bray pg 33
Cm_a    =   -2.12;               % [/rad]      Bray pg 33
Cm_adot =  -11.0;                % [/rad]      Bray pg 33
Cm_q    =  -36.6;                % [/rad]      Bray pg 33
Cm_de   =   -1.76;               % [/rad]      Bray pg 33
%Cmfade Cmfade.dat 0 1 1         %             Bray pg 34, Figure 4.2

CY_beta =   -0.819;              % [/rad]      Bray pg 33
%CY_p      <CY_p>                % [/rad]      no data
%CY_r      <CY_r>                % [/rad]      no data
%CY_da     <CY_da>               % [/rad]      no data
CY_dr   =    0.191;              % [/rad]      Bray pg 33
%CYfada <CYfada.dat>             % []          no data
%CYfbetadr CYfbetadr.dat 0 1 1   % []          Bray pg 62, Figure 4.19

Cl_beta =   -0.023;              % [/rad]      Bray pg 33
Cl_p    =   -0.450;              % [/rad]      Bray pg 33
Cl_r    =    0.265;              % [/rad]      Bray pg 33
Cl_da   =   -0.161;              % [/rad]      Bray pg 33 (sign reversed)
Cl_dr   =   -0.00229;            % [/rad]      Bray pg 33
%Clfada Clfada.dat 0 1 1         %             Bray pg 58, Table 4.8

Cn_beta =    0.109;              % [/rad]      Bray pg 33
Cn_p    =   -0.110;              % [/rad]      Bray pg 33
Cn_r    =   -0.200;              % [/rad]      Bray pg 33
Cn_da   =    0.0200;             % [/rad]      Bray pg 33 (sign reversed)
Cn_dr   =   -0.0917;             % [/rad]      Bray pg 33
%Cnfada Cnfada.dat 0 1 1         %             Bray pg 61, Table 4.9
%Cnfbetadr Cnfbetadr.dat 0 1 1   %             Bray pg 63, Figure 4.20


% Calculate CL, CD, CM

CL = CL0 + CL_a*alpha + CL_de*elevator +...
     CL_adot*alpha_dot*cbar/(2*Vt)+CL_q*q*cbar/(2*Vt);
 
CM = Cm0 + Cm_a*alpha + Cm_de*elevator +... 
     Cm_adot*alpha_dot*cbar/(2*Vt)+Cm_q*q*cbar/(2*Vt);
 
CD = CD0 + CD_a*alpha + CD_de*elevator;

% Calculate CY, Cl, CN

CY = CY_beta*beta + CY_dr*rudder;

Cl = Cl_beta*beta + Cl_da*aileron + Cl_dr*rudder + ...
    Cl_p*p*bw/(2*Vt)+Cl_r*r*bw/(2*Vt);

CN = Cn_beta*beta + Cn_da*aileron + Cn_dr*rudder + ...
    Cn_p*p*bw/(2*Vt)+Cn_r*r*bw/(2*Vt);

out(1) = CD;
out(2) = CY;
out(3) = CL;
out(4) = Cl;
out(5) = CM;
out(6) = CN;

 