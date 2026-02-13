%% 154A Project 
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

Crth = (2*Sth)/(bw*(1+lambdath)) ; % Cr at wing 
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
epsilon=0.25; % down wash 
display(hact) 
display(Lh)


% Neautral point WRT to LE of wing 
hn= (hacw + hact*((Sth/Sw)*(at/aw))*(1-epsilon))/((1)+((Sth/Sw)*(at/aw)*(1-epsilon))) ; 
display(hn) % non dimensionalized 

N= hn*Cw ; % dimensionalized neutral point 
display(N) 
% for the moment I will assume a small tail incedednce angle so it is
% negliglbe 