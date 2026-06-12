%% % MAE 171A
% MATLAB Project
% Adela Reyes 

clear, clc;

% 3.

% find poles (roots of denominator)
% denominator = (800000)s^4 + (38487000)s^3 + (1480850000)s^2 + (1375000000)s + 4e10;

poles = roots([800000 38487000 1480850000 1375000000 4e10]);
disp(poles);

% 4. 
m1 = 250;     
m2 = 50;      
b1 = 1000;    
b2 = 1500;    
k1 = 15000;   
k2 = 20000;   

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

% Step Responses
figure;
subplot(2,1,1)
step(G1)
title('Step Response: G1 (F to y)')
grid on

subplot(2,1,2)
step(G2)
title('Step Response: G2 (w to y)')
grid on

% Impulse Responses 
figure;
subplot(2,1,1)
impulse(G1)
title('Impulse Response: G1 (F to y)')
grid on

subplot(2,1,2)
impulse(G2)
title('Impulse Response: G2 (w to y)')
grid on