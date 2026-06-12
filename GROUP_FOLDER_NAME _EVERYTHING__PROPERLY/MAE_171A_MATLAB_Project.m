%% MAE 171A
% MATLAB Project
% Adela Reyes 

clear, clc, close all;

%% 3. Poles of transfer functions

m1 = 2500;     
m2 = 320;      
b1 = 350;    
b2 = 15000;    
k1 = 80000;   
k2 = 500000;   

% Denominator coefficients
den = [m1*m2, ...
       (b1*(m1+m2) + b2*m1), ...
       (b1*b2 + k1*(m1+m2) + k2*m1), ...
       (b1*k2 + b2*k1), ...
       k1*k2];

% Numerator coefficients
num1 = [(m1+m2), b2, k2];
num2 = [-m1*b2, -m1*k2, 0, 0];

% Transfer Functions
G1 = tf(num1, den);
G2 = tf(num2, den);

% Find poles
poles_G1 = pole(G1);
poles_G2 = pole(G2);

disp('Poles of G1:')
disp(poles_G1)

disp('Poles of G2:')
disp(poles_G2)

%% 4. Step Responses

figure;
subplot(2,1,1)
step(G1)
title('Step Response: G1 (F to y)')
grid on

subplot(2,1,2)
step(G2)
title('Step Response: G2 (w to y)')
grid on

%% 4. Impulse Responses 

figure;
subplot(2,1,1)
impulse(G1)
title('Impulse Response: G1 (F to y)')
grid on

subplot(2,1,2)
impulse(G2)
title('Impulse Response: G2 (w to y)')
grid on

%% 7. Bode and Nyquist plots for G1(s)

figure('Name', 'G1(s): Input F to Output y');

subplot(1,2,1);
bode(G1);
title('Bode Plot of G1(s)');
grid on;

subplot(1,2,2);
nyquist(G1);
title('Nyquist Diagram of G1(s)');
grid on;

%% 8. Bode and Nyquist plots for G2(s)

figure('Name', 'G2(s): Disturbance w to Output y');

subplot(1,2,1);
bode(G2);
title('Bode Plot of G2(s)');
grid on;

subplot(1,2,2);
nyquist(G2);
title('Nyquist Diagram of G2(s)');
grid on;

%% 9. Root Locus for G1(s)

figure;
rlocus(G1)
grid on
title('Root Locus of G1(s): Input F to Output y')
xlabel('Real Axis')
ylabel('Imaginary Axis')
axis([-30 5 -50 50])

%% Homework 7 - Problem 3

% Part (a) - Verify long oscillations with step response
figure;
step(G1,20)    
grid on
title('G1 Step Response (Long Time Horizon)')

% Display step response characteristics
info = stepinfo(G1);
disp('Step Response Characteristics:')
disp(info)

% Part (b) - Design specifications region
OS = 5;         
Ts = 1;          

% Damping ratio from overshoot
zeta = -log(OS/100)/sqrt(pi^2 + (log(OS/100))^2);

% For 1% settling time
sigma = 4.6/Ts;

disp(['Required damping ratio zeta = ', num2str(zeta)])
disp(['Required real part <= -', num2str(sigma)])

figure;
rlocus(G1)
grid on
hold on

% Overlay damping ratio lines
sgrid(zeta,[])

% Settling time boundary
xline(-sigma,'r--','T_s = 1 s')

title('Root Locus with 5% OS and 1 s Settling Time Constraints')
xlabel('Real Axis')
ylabel('Imaginary Axis')
axis([-30 5 -50 50])

%% Homework 8 
% Part 1

C = 1;   % unity controller

% Reference to output
T_r2y = feedback(C*G1,1);

figure
step(T_r2y)
grid on
title('Step Response: Reference to Output')