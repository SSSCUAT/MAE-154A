

Wguess = 1000;
W = Wguess; 

for i = 1:10
Wto(i) = W;



%S = 7.12;     % 

S = 18; %W/(.5*.00238*60^2 * 1.2 );

A = 12; %8.68;
Sh = 3.63;    
Sv = 1.51;
bh = 3.88;
bv = 1.62;


hac = 0.25;
c = 0.91;
ch = 0.94;
cv = .93;
hach = 0.25;
hacv = 0.25;


%% Component Weight Estimates- Nicolai

%% Wing Weight

N=6.6;              %Ultimate Load Factor (1.5 times limit load factor)(GIVEN)
Delta=0*pi/180;%Deg %Wing 1/4 chord sweep angle
tr=1;               %Taper Ratio
tc=0.12;            %Maximum Thickness Ratio (GIVEN)
Ve=80;%kts         %Equivalent Vmax at SL

Ww=96.948*((W*N/10^5)^0.65*(A/cos(Delta))^0.57*(S/100)^0.61*((1+tr)/(2*tc))^0.36*(1+Ve/500)^0.5)^0.993

%% Fuselage Weight

lf=7.333; %ft       %Fuselage Length
WF=1.75; %ft        %Fuselage Width
D=1/12; %ft          %Fuselage Max Depth

Wf=200*((W*N/10^5)^0.286*(lf/10)^0.857*((WF+D)/10)*(Ve/100)^0.338)^1.1

%% Horizontal Tail Weight



lh=35 / 12 + (.5 - hac) * c - (.5 - hach) * ch; %ft       %Distance from Wing MAC to Tail MAC
thr=ch*.12*12; %inches      %horizontal tail max root thickness (chord * thick/chord)

Wht=127*((W*N/10^5)^0.87*(Sh/100)^1.2*(lh/10)^0.483*(bh/thr)^0.5)^0.458

%% Vertical Tail Weight

tvr=cv*.12*12; %in    %Vertical Tail Max Root Thickness (chord * thick/chord * in/ft)

Wvt= (2)*  98.5*((W*N/10^5)^0.87*(  (.5)*  Sv/100)^1.2*(  (.5)*  bv/tvr)^0.5)^0.458

%% Landing Gear Weight

%Llg=18; %in     %Length of Main Landing Gear Strut
%Nland=2;        %Ultimate Load Factor at Wland
%Wlg=0.054*(Llg)^0.501*(W*Nland)^0.684

%don't need niccolai if we have specific landing gear 
Wlg = 1.5;

%% TOTAL STRUCTURAL WEIGHT

Wstruct=Ww+Wf+Wht+Wvt+Wlg


%% Total Propulsion Unit (minus Fuel system) Weight

Weng=3.3; %(lbs)     %Bare Engine Weight
Neng=1;             %# Engines

Wp=2.575*(Weng)^0.922*Neng    %this equation likely over-estimates propulsion unit weight for a small UAV


%% Fuel Weight

Wfu = 5;   %(lbs)  

%% Fuel System Weight

%rhof = 6.739; %lb/gal fuel mass density JP-8
%Fg = Wfu / rhof; %gal               %Total Fuel
%tankint=1; %percent         %Percent of Fuel Tanks that are integral
%Nt=2;                         %Number of Separate Fuel Tanks
%Wfs=2.49*((Fg)^0.6*(1/(1+tankint))^0.3*Nt^0.2*Neng^0.13)^1.21

% specific fuel system weights (fuel tanks, lines) likely can be found for your aircraft, if so, use those actual values instead of the niccolai equations.
Wfs = .25;

%% Surface Controls Weight

Wsc=1.066*W^0.626;  

%% Avionics Weight - use weights of specific sensors you choose

Wau=1       


%% Payload Weight

Wpl=20;

%% TOTAL WEIGHT

Wto(i)=Wstruct+Wp+Wfs+Wsc+Wpl+Wfu
W = Wto(i)

end
%figure; grid on;
hold on

plot([Wguess Wto],'.-m')
