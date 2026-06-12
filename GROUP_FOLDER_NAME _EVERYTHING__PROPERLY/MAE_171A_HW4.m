%% MAE 171 
% Homework 4 
% Adela Reyes


%% Problem 1

clear, clc ;
s = tf('s');

% a)
G1a = 1/(s+10);

figure;
nyquist(G1a)
grid on
title('Problem 1(a): G(s) = 1/(s+10)')

% b)
G1b = (s+1)/((s+3)*(s+10));

figure;
nyquist(G1b)
grid on
title('Problem 1(b): G(s) = (s+1)/((s+3)(s+10))')

% c)
G1c = (1/s)*exp(-2*s);

figure;
nyquist(G1c)
grid on
title('Problem 1(c): G(s) = (1/s)e^{-2s}')

%% Problem 2

% a)
G2a = (s+2)/(s+10);

figure;
nyquist(G2a)
grid on
title('Problem 2(a): G(s) = (s+2)/(s+10)')

% b) 
G2b = 1/((s+10)*(s+2)^2);

figure;
nyquist(G2b)
grid on
title('Problem 2(b): G(s) = 1/((s+10)(s+2)^2)')


% c)
G2c = ((s+10)*(s+1))/((s+100)*(s+2)^3);

figure;
nyquist(G2c)
grid on
title('Problem 2(c): G(s) = ((s+10)(s+1))/((s+100)(s+2)^3)')