%UAV sim -- D.Toohey

function [sys,x0,str,ts] = aero_dyn(t,x,u,flag)

switch flag,
  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  % Initialize the states, sample times, and state ordering strings.
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=mdlDerivatives(t,x,u);
    
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  % Return the outputs of the S-function block.
  case 3
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  % There are no termination tasks (flag=9) to be handled.
  % Also, there are no continuous or discrete states,
  % so flags 1,2, and 4 are not used, so return an emptyu
  % matrix 
  case { 2, 4, 9 }
    sys=[];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unexpected flags (error handling)%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Return an error message for unhandled flag values.
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end


%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array
%
sizes = simsizes;

sizes.NumContStates  = 13;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 14;
sizes.NumInputs      = 12;
sizes.DirFeedthrough = 0;     % 1 : yes, 0 : no
sizes.NumSampleTimes = 1;     % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%

init_dynamics;

Phi = Eulini(1);
Theta = Eulini(2);
Psi = Eulini(3);
   
q0i=cos(Phi/2)*cos(Theta/2)*cos(Psi/2) + sin(Phi/2)*sin(Theta/2)*sin(Psi/2);
q1i=sin(Phi/2)*cos(Theta/2)*cos(Psi/2) - cos(Phi/2)*sin(Theta/2)*sin(Psi/2);
q2i=cos(Phi/2)*sin(Theta/2)*cos(Psi/2) + sin(Phi/2)*cos(Theta/2)*sin(Psi/2);
q3i=cos(Phi/2)*cos(Theta/2)*sin(Psi/2) - sin(Phi/2)*sin(Theta/2)*cos(Psi/2);

% Initial State vector
x0=[norm(Velini) alpha0 beta0  q0i q1i q2i q3i Angvelini Posini]';

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];      % continuous sample time

% end mdlInitializeSizes


%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

%init_dynamics;

weight = 420; 
grav = 32.17405;
mass = weight/grav; % [lbm]
I_xx = 34.832;    % [slug-ft^2] IAI, MSS notes
I_yy = 67.08;     % [slug-ft^2] IAI, MSS notes
I_zz = 82.22;     % [slug-ft^2] IAI, MSS notes
I_xz = -4.902;    % [slug-ft^2] IAI, MSS notes
Inertia = [I_xx 0 I_xz; 0 I_yy 0; I_xz 0 I_zz];


% aerodynamic Forces and Moments from aero
Txb = u(1);
Tyb = u(2);
Tzb = u(3);
Fxs = u(4);  % Aerodynamic Forces
Fys = u(5);
Fzs = u(6);
Ls = u(7);   % Aerodynamic Moments
Ms = u(8);
Ns = u(9);
Gxe = u(10); % Gravity Accel. in Earth Frame (North, East, Down)
Gye = u(11);
Gze = u(12);

Ta = [Txb Tyb Tzb]'./mass;   % Thrust accel in the body frame
G = [Gxe Gye Gze]';

%state variables
Vt=x(1);
alpha=x(2);
beta=x(3);
q0=x(4);
q1=x(5);
q2=x(6);
q3=x(7);
P=x(8);
Q=x(9);
R=x(10);
pN=x(11);
pE=x(12);
pD=x(13);


%check normalization of quaternion
q_mag = (q0^2+q1^2+q2^2+q3^2)^.5;

q0 = q0/q_mag;
q1 = q1/q_mag;
q2 = q2/q_mag;
q3 = q3/q_mag;


%Matrix for crossproduct with Angular Velocity of a body
OMEGA3=[
    0 -R  Q;
    R  0 -P;
    -Q  P  0];

%Matrix for Kinematic Relation when using Quaternions
OMEGA4=[
    0  P  Q  R;
    -P  0 -R  Q;
    -Q  R  0 -P
    -R -Q  P  0];

%Rotational Transformation Matrix from Body to Earth
C11=q0^2+q1^2-q2^2-q3^2;  C12=2*(q1*q2-q0*q3);      C13=2*(q1*q3+q0*q2);
C21=2*(q1*q2+q0*q3);      C22=q0^2-q1^2+q2^2-q3^2;  C23=2*(q2*q3-q0*q1);
C31=2*(q1*q3-q0*q2);      C32=2*(q2*q3+q0*q1);      C33=q0^2-q1^2-q2^2+q3^2;

C_b2e=[
    C11 C12 C13;
    C21 C22 C23;
    C31 C32 C33];

%Conversion of aerodynamic forces and moments from stability axis to body axis
Fx=Fxs*cos(alpha)-Fzs*sin(alpha);
Fy=Fys;
Fz=Fxs*sin(alpha)+Fzs*cos(alpha);
L=Ls*cos(alpha)-Ns*sin(alpha);
M=Ms;
N=Ls*sin(alpha)+Ns*cos(alpha);
%Conversion:Vt,alpha,beta --> U,V,W
U=Vt*cos(alpha)*cos(beta);
V=Vt*sin(beta);
W=Vt*sin(alpha)*cos(beta);
vB=[U V W]';
%Perfect Accelerometer Outputs
A_x=Fx/mass;
A_y=Fy/mass;
A_z=Fz/mass;
A=[A_x A_y A_z]';

%Force Equation in body axis
temp=A + Ta + C_b2e'*G - OMEGA3*vB;
UDOT=temp(1);
VDOT=temp(2);
WDOT=temp(3);

%Conversion:U,V,W-->Vt,alpha,beta
sys(1)=(U*UDOT+V*VDOT+W*WDOT)/Vt;
sys(2)=(U*WDOT-W*UDOT)/(U*U+W*W);
sys(3)=(Vt*VDOT-V*(sys(1)))/cos(beta)/Vt^2;

%Kinematic Equation
temp=-0.5*OMEGA4*[q0 q1 q2 q3]';
sys(4)=temp(1);
sys(5)=temp(2);
sys(6)=temp(3);
sys(7)=temp(4);

%Moment Equation in body axis
Lpp = L + Inertia(1,3)*P*Q - (Inertia(3,3)-Inertia(2,2))*R*Q;
Np = N - (Inertia(2,2)-Inertia(1,1))*P*Q - Inertia(1,3)*R*Q;


sys(8)= (Lpp*Inertia(3,3) - Np*Inertia(1,3))/(Inertia(1,1)*Inertia(3,3) - Inertia(1,3)^2);
sys(9)= (M - (Inertia(1,1) - Inertia(3,3))*P*R - Inertia(1,3)*(P^2 - R^2))/Inertia(2,2);
sys(10)= (Np*Inertia(1,1) + Lpp*Inertia(1,3))/(Inertia(1,1)*Inertia(3,3) - Inertia(1,3)^2);

%Navigation Equation
temp=C_b2e*vB;
sys(11)=temp(1);
sys(12)=temp(2);
sys(13)=temp(3);

%sys         (1) (2)  (3)  (4) (5) (6) (7) (8)(9)(10)(11)(12) (13)
%derivate of Vt alpha beta  q0  q1  q2  q3  P  Q  R  pN   pE   pD

% end mdlDerivatives



%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)


if x(13) > 0
   x(13) = 0; 
end

    sys=[x', t];
   
% end mdlOutputs