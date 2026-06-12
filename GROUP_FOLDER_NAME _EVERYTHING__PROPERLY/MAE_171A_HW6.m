%% MAE 171A Homework 6
% Adela Reyes

clear, clc;

%% Problem 2 (5.5)

s = tf('s');

% a) 
La = 1/(s^2 + 3*s + 10);

figure;
rlocus(La);
grid on;
title('Root Locus for Part (a)');

% b) 
Lb = 1/(s*(s^2 + 3*s + 10));

figure;
rlocus(Lb);
grid on;
title('Root Locus for Part (b)');

% c)
Lc = (s^2 + 2*s + 8)/(s*(s^2 + 2*s + 10));

figure;
rlocus(Lc);
grid on;
title('Root Locus for Part (c)');

% d)
Ld = (s^2 + 2*s + 12)/(s*(s^2 + 2*s + 10));

figure;
rlocus(Ld);
grid on;
title('Root Locus for Part (d)');

%% Problem 3 (5.7)

clear; clc;

s = tf('s');

% Part (b)
Lb = (s+3)/(s^2*(s+10)*(s^2+6*s+25));

figure;
rlocus(Lb);
grid on;
axis([-15 10 -15 15]);
title('Root Locus for 5.7(b)');


% Part (c)
Lc = (s+3)^2/(s^2*(s+10)*(s^2+6*s+25));

figure;
rlocus(Lc);
grid on;
axis([-15 10 -15 15]);
title('Root Locus for 5.7(c)');


% Part (e)
Le = ((s+1)^2 + 1)/(s^2*(s+2)*(s+3));

figure;
rlocus(Le);
grid on;
axis([-15 10 -15 15]);
title('Root Locus for 5.7(e)');