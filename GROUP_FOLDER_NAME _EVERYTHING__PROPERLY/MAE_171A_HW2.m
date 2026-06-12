% MAE 171A
% Homework 2
% Adela Reyes 

clear, clc;

%% Problem 4

% a)
% find poles (roots of denominator)
% denominator = s^3 + s^2 + 3*s + 2;

poles = roots([1 1 3 2]);
disp(poles);

% b) 
% denominator of c/r = s^3 + s^2 + 3*s + 3;

poles = roots([1 1 3 3]);
disp(poles);

% c) 
% denominator of c/r = s^3 + s^2 + 3*s + 7 , where Gc = 5;

poles = roots([1 1 3 7]);
disp(poles);
