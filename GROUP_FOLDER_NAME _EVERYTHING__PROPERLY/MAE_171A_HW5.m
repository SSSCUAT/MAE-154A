%% MAE 171A Homework 5

% Problem 3 
% b)

clc, clear;

s = tf('s');

% Parameters
K = 13.3;
tau1 = 1;
taum = 0.5;
J = 4;
b = 1;

% Loop transfer function
L = K / ((tau1*s + 1)*(J*s + b)*(taum*s + 1));

% Gain Margin and Phase Margin
figure;
margin(L)
grid on
title('Gain Margin and Phase Margin')