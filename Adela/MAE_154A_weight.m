clear, clc

Wguess = 20;
W = Wguess;
for i = 1:10
    Wto(i) = W;

    Arw = 7.4;         % aspect ratio]
    CL_MAX = 1.5;
    rho = 0.0023769;
    V_stall = 55;
    Sw = (W*2) / (CL_MAX*rho*(V_stall^2)); % wing area
    bw = sqrt(Arw * Sw); % wing span (tip to tip)
    Sth = 0.2 *Sw;
    bth = 1.5;
    S_vt = .41;
    b_vt = 0.7;
    c_vt = .5;

    c_ht = 0.38;
    hac_ht = 0.25;
    hac_vt = 0.25;
    c = 0.8;
    hac = 0.25;
    
% Component Weight Estimates- Nicolai
% Wing Weight 
    
    N=6.6; %Ultimate Load Factor (1.5 times limit load factor)(GIVEN)
    Delta=0*pi/180;%Deg %Wing 1/4 chord sweep angle
    lambdaw =0.4; %Taper Ratio
    tc=0.12; %Maximum Thickness Ratio (GIVEN)
    V_max_kts=80;%kts %Equivalent Vmax at SL

    Ww = 96.948*((W*N/10^5)^0.65*(Arw/cos(Delta))^0.57*(Sw/100)^0.61*((1+lambdaw)/(2*tc))^0.36*(1+V_max_kts/500)^0.5)^0.993;

%% Fuselage Weight

     L_fuse = 5; %ft %Fuselage Length
     W_fuse =.67; %ft %Fuselage Width
     D_fuse = 1/24; %ft %Fuselage Max Depth

    Wf=200*((W*N/10^5)^0.286*(L_fuse/10)^0.857*((W_fuse+D_fuse)/10)*(V_max_kts/100)^0.338)^1.1;

%% Horizontal Tail Weight

    lh=35 / 12 + (.5 - hac) * c - (.5 - hac_ht) * c_ht; %ft %Distance from Wing MACto Tail MAC
    thr=c_ht*.12*12; %inches %horizontal tail max root thickness (chord * thick/chord)

    Wht=127*((W*N/10^5)^0.87*(Sth/100)^1.2*(lh/10)^0.483*(bth/thr)^0.5)^0.458;

%% Landing Gear Weight

% Llg=18; %in    %Length of Main Landing Gear Strut
% Nland=2; %Ultimate Load Factor at Wland
% Wlg=0.054*(Llg)^0.501*(W*Nland)^0.684

%don't need niccolai if we have specific landing gear
    Wlg = 1.5;
%% Weight of propellar 
    Wprop = 0.24;

%% Vertical Tail Weight

    tvr=c_vt*.12*12; %in %Vertical Tail Max Root Thickness (chord * thick/chord * in/ft)

    Wvt= (2)* 98.5*((W*N/10^5)^0.87*( (.5)* S_vt/100)^1.2*( (.5)* b_vt/tvr)^0.5)^0.458;

%% TOTAL STRUCTURAL WEIGHT
    Wstruct=Ww+Wf+Wht+Wvt+Wlg;

%% Total Propulsion Unit (minus Fuel system) Weight

    Weng=1.76; %(lbs) %Bare Engine Weight
    Neng=1; %# Engines

    Wp=2.575*(Weng)^0.922*Neng; %this equation likely over-estimates propulsion unit weight for a small UAV

%% Fuel Weight

     Wfuel = 1.5; %(lbs)

%% Fuel System Weight

%rhof = 6.739; %lb/gal fuel mass density JP-8
%Fg = Wfu / rhof; %gal %Total Fuel
%tankint=1; %percent %Percent of Fuel Tanks that are integral
%Nt=2; %Number of Separate Fuel Tanks
%Wfs=2.49*((Fg)^0.6*(1/(1+tankint))^0.3*Nt^0.2*Neng^0.13)^1.21

% specific fuel system weights (fuel tanks, lines) likely can be found for your aircraft, if so, use those actual values instead of the niccolai equations.
    Wfs = .25;
%% Surface Controls Weight

    Wsc=1.066*W^0.626;
%% Avionics Weight - use weights of specific sensors you choose (camera, servos, computer, GPS & Battery)
     Wau=1.62;
%% Payload Weight
     Wpay=2;
%% TOTAL WEIGHT
    Wto(i)=Wstruct+Wp+Wfs+Wsc+Wpay+Wfuel;
    W = Wto(i);

    Ww_arr(i) = Ww;
    Wf_arr(i) = Wf;
    Wht_arr(i) = Wht;
    Wvt_arr(i) = Wvt;
    Wstruct_arr(i) = Wstruct;
    Wp_arr(i) = Wp;
    Wau_arr(i) = Wau;
end

% W at elbow!!
index = 3;            % should choose x=4
W_pick = Wto(index)

% Outputs of component weights @ index!!
Ww = Ww_arr(index)
Wf = Wf_arr(index)
Wht = Wht_arr(index)
Wvt = Wvt_arr(index)
Wstruct = Wstruct_arr(index)
Wp = Wp_arr(index);
Wau = Wau_arr(index);

figure; grid on; hold on
plot([Wguess Wto], '.-m')
plot(index+1, W_pick, 'or', 'MarkerSize', 10, 'LineWidth', 2)  % +1 bc Wguess is the first point

