%% 154A Project 
%known values 

% Aerodynamic Variable Definitions
CL = 1.5;    % 3D wing lift coefficient
AR = 7.4;    % Wing aspect ratio (Example value for Cessna 182)
CD_p = 0.027; %parasite Drag 
e = 0.80;            % Oswald's efficiency factor

% Induced Drag Calculation
% Formula from image: CDi = CL^2 / (pi * AR)
CDi = (CL^2) / (pi * AR);

% Display Result
fprintf('The Induced Drag Coefficient (CDi) is: %.4f\n', CDi);


%L=W=0.5*rho*v^2*S*CL 
W=20;   % Estimate lbs
V=55;    %ft/s
b=5.67 ; 
S=4.35 ; 
rho= 0.002378;  %slugs/ft³
Cl= W./(0.5*rho*S*V.^2) ; 
disp(Cl)
%% 
clear
clc
%all units should be in ft 
Arw=7.4 ; 
%bw = 4.5; 
bw =5.617 ; 
%Sw= (bw.^2)/Arw; 
Sw = 4.264931025 ; 
Sth = 0.5331163781 ; %surface area of horizontal tail 
bth = 1.460296378 ; 


%size of wings and horizontal tail 
lambdaw= 0.35 ; % taper ratio of wing 
lambdat= 0.5 ; % taper ratio of tail (Horizontal) 
%Sw= bw/2 *( Cr+ Ctip) %Cr=chord at root, bw=span %Ctip=LambaCr 
Crw = (2*Sw)/(bw*(1+lambdaw)) ; % Cr at wing 
Ctw= lambdaw* Crw ; % chord length at tip 
Cw = (2/3) * Crw *((1 + lambdaw + lambdaw.^2)/(1+lambdaw)) ; %effective chord length
display(Crw)
display(Ctw)
display(Cw) 




%% 
%Static Margin 



Cw= 0.867; %chord of wing 
Ct=0.40  ; %chord of tail 
e= 0.80 ; % oswald efficieny 
% im assuming a 2d lifting slope of 2pi from TAT 
cla= 2*pi ; % 2d lift curve slope of NACA 2412 
Arw= 7.4 ; 
Art= 4 ; 

at = cla/ (1+ (cla/(pi*Art*e))) ; 
aw = cla/ (1+ (cla/(pi*Arw*e))) ; 
display(at) 
display(aw)  



hacw= 0.25 ; % will always be at 0.25 quarter chord of wing 
hact= 4.29 ; % distance from LE of wing to AC of Tail 
St=0.695 ; % surface area of tail 
Sw= 5.57 ; % surcae area of wing 
epsilon=0.25; % down wash 


hn= (hacw + hact*((St/Sw)*(at/aw))*(1-epsilon))/((1)+((St/Sw)*(at/aw)*(1-epsilon))) ; 
display(hn) 
