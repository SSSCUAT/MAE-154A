%% 154A Project 
%(do not refrence this code it is for me) 
clear
clc
%all units should be in ft 
Arw=7.4 ; 
bw =5.617 ; 
%Sw= (bw.^2)/Arw; 

Sw = 4.264931025 ; % surface area of wing 
Sth = 0.5331163781 ; %surface area of horizontal tail 
bth = 1.460296378 ;  % span of horizontal tail 


%size of wings and horizontal tail 
lambdaw= 0.40 ; % taper ratio of wing 
lambdath= 0.5 ; % taper ratio of tail (Horizontal) 
%Sw= bw/2 *( Cr+ Ctip) %Cr=chord at root, bw=span %Ctip=LambaCr 
Crw = (2*Sw)/(bw*(1+lambdaw)) ; % Cr at wing 
Ctw= lambdaw* Crw ; % chord length at tip 
Cw = (2/3) * Crw *((1 + lambdaw + lambdaw.^2)/(1+lambdaw)) ; %effective chord length
display(Crw)
display(Ctw)
display(Cw) 

Crth = (2*Sth)/(bth*(1+lambdath)) ; % Cr at wing 
Ctth= lambdath* Crth; % chord length at tip 
Cth = (2/3) * Crth *((1 + lambdath + lambdath.^2)/(1+lambdath)) ; %effective chord length
display(Crth)
display(Ctth)
display(Cth)




%Static Margin 

e= 0.85 ; % oswald efficieny (guess) 


cla= 5.73  ; % 2d lift curve slope of NACA 2412 units= rad^-1 
% data from virigna tech unversity 


Art= 4 ; % static margin of tail 

at = cla / (1+ (cla/(pi*Art*e))) ; % 3d lifitng slope of tail
aw = cla/ (1+ (cla/(pi*Arw*e))) ; % 3d lifitng slope of wing 
display(at) 
display(aw)  



hacw= 0.25 ; % will always be at 0.25 quarter chord of wing 
Lh= 4*Cw ; %distance from wing aerodynamic center (about quarter chord) 
% to horizontal tail aerodynamic center

xacw = hacw*Cw;  % quater chord of wing in ft 
xact= xacw+ Lh ; % distance of ac of tail wrt to LE of wing 
hact= xact/Cw ; % distance from LE of wing to AC of Tail 
downwash=0.25; % down wash 
display(hact) 
display(Lh)


% Neautral point WRT to LE of wing 
hn= (hacw + hact*((Sth/Sw)*(at/aw))*(1-downwash))/((1)+((Sth/Sw)*(at/aw)*(1-downwash))) ; 
display(hn) % non dimensionalized 

N= hn*Cw ; % dimensionalized neutral point 
display(N) 
% for the moment I will assume a small tail incedednce angle so it is
% negliglbe or tail incidence (it) ~-1degree if necessary for now 
it= -1 * pi/180 ; % incidence angle equals -1degree but I calculated it in rads
%hcg=? ; %distance from leading edge of wing to cg normalized by dividing by wing span (Cw) 

% I am going to continue assuming a downwash of 0.25 from cessna data 
Cmacw = -0.05;  %from data accroding to virnigna tech libraries 
Cmcgw= Cmacw + aw*alpha(hcg-hacw); %can be calcualted function of (alpha) 
   % for Cmcg,t we need lt lt= distance from cg to tail ac 
Cmcgt= -( VH* at*(1-downwash) )*alpha + at*VH*it; % (function of alpha) 
VH= (lt/c) * (St/Sw); %tail volume ratio 
Cmcg = Cmcgw + Cmcgt ; 
Clalpha = aw + at*(Sth/ Sw)* (1-downwash) ; 
Cmalpha= -CLalpha (hn-hcg) ; 
% tau= elevator effectivness for now I will guess an effectivness of 0.5 
tau= 0.5 ; 
% deltae= elevator deflection 
Cm0= Cmacw + VH*at*it ; 
Cl0= -at * (Sth/SW) * it ; 
Cldeltae= at* (Sth/Sw) * tau ;   %( not sure if positive or negative ) pretty sure its positive 
Cmdeltae= - at*VH * tau ; %( not sure if positive or negative ) pretty sure its negative 
deltae= - ((Cm0*Clalpha) + (Cmalpha*Cl) )./ ((Clalpha*Cmdeltae - Cmalpha * Cldeltae)) ; 
%deltae is a function of Cl 


%% MAE 154A - Longitudinal Stability / Trim Setup (TEAM REFERENCE)
% -------------------------------------------------------------------------
% PURPOSE:
%   1) Compute wing & tail geometry (root/tip/MAC)
%   2) Compute 3D lift-curve slopes (aw, at)
%   3) Compute neutral point (hn) and static margin (SM)
%   4) Set up the coefficients needed for trim / elevator deflection
%
%  CONVENTION:
%   All "h" quantities are measured from the WING LE and non-dimensionalized
%   by the WING MAC (Cw).
% Units:
%   length = ft
%   angles = rad
%   aerodynamic slopes = 1/rad
%   coefficients (CL, Cm, etc.) are nondimensional
% -------------------------------------------------------------------------

clear; clc;

%
% 1) INPUTS YOU SHOULD EDIT (KNOWN VALUES / ASSUMPTIONS)

% -----------------------------
% Wing geometry (KNOWN)
% -----------------------------
Arw     = 7.4;             % Wing aspect ratio [-]
bw      = 5.617;           % Wing span [ft]
Sw      = 4.264931025;     % Wing planform area [ft^2]
lambdaw = 0.40;            % Wing taper ratio ct/cr [-]

% -----------------------------
% Horizontal tail geometry (KNOWN)
% -----------------------------
Sth      = 0.5331163781;   % Horizontal tail area [ft^2]
bth      = 1.460296378;    % Horizontal tail span [ft]
lambdath = 0.50;           % Tail taper ratio ct/cr [-]
Art      = 4;              % Tail aspect ratio [-] (given/assumed)

% -----------------------------
% Aero assumptions (EDIT AS NEEDED)
% -----------------------------
e        = 0.85;           % Oswald efficiency factor [-] (guess)
cla      = 5.73;           % 2D airfoil lift curve slope [1/rad] (check units!)
downwash  = 0.25;           % dε/dα downwash gradient [-] (assumed)
it       = -1*pi/180;      % Tail incidence angle [rad] (-1 deg)
Cmacw    = -0.05;          % Wing Cm about wing AC (approx from airfoil data)
tau      = 0.5;            % Elevator effectiveness τe [-] (guess)

% -----------------------------
% Layout / reference positions (KNOWN FROM CG BUILD-UP)
% -----------------------------
x_wle    = 0.75;           % Wing LE location from nose [ft] (given)
x_cg_str = 1.7147;         % Structural CG from nose [ft]
x_cg_tot = 1.4103;         % CG with avionics from nose [ft]
x_cg_0   = 1.5071;         % Initial CG w/ payload+fuel from nose [ft]

% -----------------------------
% Tail placement ASSUMPTION
% -----------------------------
% Lh is distance from wing AC to tail AC (this is a design/layout choice)
% You currently assume Lh = 4*Cw (updated after Cw is computed).
% -------------------------------------------------------------------------


% ========================================================================
% 2) PLANFORM GEOMETRY (ROOT, TIP, MAC) - WING + TAIL
%

% ---- Wing chords ----
Crw = (2*Sw)/(bw*(1 + lambdaw));  % Wing root chord [ft]
Ctw = lambdaw * Crw;              % Wing tip chord [ft]

% Wing MAC (mean aerodynamic chord) for trapezoid
Cw  = (2/3) * Crw * ((1 + lambdaw + lambdaw^2)/(1 + lambdaw));  % [ft]

fprintf("\n--- WING GEOMETRY ---\n");
fprintf("Crw (root chord)  = %.4f ft\n", Crw);
fprintf("Ctw (tip chord)   = %.4f ft\n", Ctw);
fprintf("Cw  (MAC)         = %.4f ft\n", Cw);

% ---- Horizontal tail chords ----
Crth = (2*Sth)/(bth*(1 + lambdath));  % Tail root chord [ft]
Ctth = lambdath * Crth;              % Tail tip chord [ft]
Cth  = (2/3) * Crth * ((1 + lambdath + lambdath^2)/(1 + lambdath)); % Tail MAC [ft]

fprintf("\n--- HORIZONTAL TAIL GEOMETRY ---\n");
fprintf("Crth (root chord) = %.4f ft\n", Crth);
fprintf("Ctth (tip chord)  = %.4f ft\n", Ctth);
fprintf("Cth  (MAC)        = %.4f ft\n", Cth);


% ========================================================================
% 3) 3D LIFT CURVE SLOPES (FINITE WING CORRECTION)

% 3D lift slope: a = a0 / (1 + a0/(pi*e*AR))
aw = cla / (1 + (cla/(pi*e*Arw)));   % Wing 3D slope [1/rad]
at = cla / (1 + (cla/(pi*e*Art)));   % Tail 3D slope [1/rad]

fprintf("\n--- LIFT CURVE SLOPES ---\n");
fprintf("aw (wing) = %.4f 1/rad\n", aw);
fprintf("at (tail) = %.4f 1/rad\n", at);


% ========================================================================
% 4) LONGITUDINAL REFERENCE LOCATIONS (MEASURED FROM WING LE)

% Wing aerodynamic center location (fraction of MAC)
hacw = 0.25;               % wing AC at 25% MAC (standard subsonic assumption)

% Wing AC measured from wing LE (dimensional)
x_ac_w_LE = hacw*Cw;        % [ft] measured from wing LE

% Tail location assumption: wing AC -> tail AC distance
Lh = 4*Cw;                  % [ft] (your assumed tail arm from wing AC)

% Tail AC measured from wing LE
x_ac_t_LE = x_ac_w_LE + Lh; % [ft] measured from wing LE

% Nondimensional tail AC location from wing LE
h_ac_t = x_ac_t_LE / Cw;    % [-]

fprintf("\n--- REFERENCE LOCATIONS (FROM WING LE) ---\n");
fprintf("Wing AC from wing LE  x_ac_w_LE = %.4f ft  (hacw = %.2f)\n", x_ac_w_LE, hacw);
fprintf("Tail AC from wing LE  x_ac_t_LE = %.4f ft\n", x_ac_t_LE);
fprintf("Tail AC nondim (h_ac_t)         = %.4f\n", h_ac_t);


% ========================================================================
% 5) NEUTRAL POINT (hn) + STATIC MARGIN (SM)

hn= (hacw + h_ac_t*((Sth/Sw)*(at/aw))*(1-downwash))/((1)+((Sth/Sw)*(at/aw)*(1-downwash))) ; 

% Compute hcg for different loading cases (nondimensional from wing LE)
% Here we convert nose-station CG -> wing-LE-based coordinate then divide by MAC
hcg_str = (x_cg_str - x_wle)/Cw;
hcg_tot = (x_cg_tot - x_wle)/Cw;
hcg_0   = (x_cg_0   - x_wle)/Cw;

% Static margin for each case
SM_str = hn - hcg_str;
SM_tot = hn - hcg_tot;
SM_0   = hn - hcg_0;

fprintf("\n--- NEUTRAL POINT & STATIC MARGIN ---\n");
fprintf("hn (neutral point) = %.4f  [-]\n", hn);
fprintf("hcg_str = %.4f,  SM_str = %.4f  (%.1f%% MAC)\n", hcg_str, SM_str, 100*SM_str);
fprintf("hcg_tot = %.4f,  SM_tot = %.4f  (%.1f%% MAC)\n", hcg_tot, SM_tot, 100*SM_tot);
fprintf("hcg_0   = %.4f,  SM_0   = %.4f  (%.1f%% MAC)\n", hcg_0,   SM_0,   100*SM_0);

% NOTE:
% If SM is negative -> CG is behind neutral point -> statically unstable.
% If SM is positive -> statically stable.

% ========================================================================
% 6) TAIL ARM lt (CG -> TAIL AC) AND TAIL VOLUME VH


% We compute tail AC from nose, then subtract CG from nose.
% This gives lt in feet (moment arm for tail forces about the CG).

% Wing AC from nose:
x_wac = x_wle + x_ac_w_LE;     % = wing LE from nose + (0.25*MAC from LE)

% Tail AC from nose:
x_tac = x_wac + Lh;            % wing AC -> tail AC

% Tail arm (CG -> tail AC) for each loading case
lt_str = x_tac - x_cg_str;
lt_tot = x_tac - x_cg_tot;
lt_0   = x_tac - x_cg_0;

% Tail volume coefficient VH = (St * lt) / (Sw * Cw)
VH_str = (Sth*lt_str)/(Sw*Cw);
VH_tot = (Sth*lt_tot)/(Sw*Cw);
VH_0   = (Sth*lt_0)  /(Sw*Cw);

fprintf("\n--- TAIL ARM & TAIL VOLUME ---\n");
fprintf("x_tac (tail AC from nose) = %.4f ft\n", x_tac);
fprintf("lt_str = %.4f ft,  VH_str = %.4f\n", lt_str, VH_str);
fprintf("lt_tot = %.4f ft,  VH_tot = %.4f\n", lt_tot, VH_tot);
fprintf("lt_0   = %.4f ft,  VH_0   = %.4f\n", lt_0,   VH_0);


% ========================================================================
% 7) LONGITUDINAL COEFFICIENTS 
% To compute trim elevator deflection you MUST choose:
%   alpha = trim angle of attack [rad]   (or use the combined trim formula with CL)
%   CL    = required lift coefficient at the trim condition
%
% These are placeholders so the script doesn't crash.
% Uncomment and set them when you're ready.

% alpha = 0;      % [rad] placeholder
% CL    = 0.7;    % [-] placeholder

% ---- Example: use structural case VH_str and hcg_str ----
% (If you want a different loading case, swap _str for _tot or _0)

% Tail volume coefficient for chosen case:
VH = VH_tot;
hcg = hcg_tot;

% Lift slope of aircraft (wing + tail contribution)
CLalpha = aw + at*(Sth/Sw)*(1 - downwash);

% Pitching moment slope using neutral point relation (lecture form)
Cmalpha = -CLalpha*(hn - hcg);

% Cm0 (from lecture note box)
Cm0 = Cmacw + VH*at*it;

% Control derivatives (first-pass)
CLdeltae = at*(Sth/Sw)*tau;   % usually positive
Cmdeltae = -at*VH*tau;        % usually negative

% Trim elevator equation (requires CL defined!)
%deltae = -((Cm0*CLalpha) + (Cmalpha*CL)) / ((CLalpha*Cmdeltae) - (Cmalpha*CLdeltae));
% fprintf("\nTrim elevator deltae = %.3f rad (%.2f deg)\n", deltae, rad2deg(deltae));


display(CLalpha)