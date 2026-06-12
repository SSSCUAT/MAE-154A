%% MAE 171A
% Homework 3
% Adela Reyes

%% Problem 4
Kvals = [-1 0 1 3 3.885 5];

figure; hold on; grid on;
xlabel('Real Axis'); ylabel('Imag Axis');
title('Roots in the s-plane for different K values');

for K = Kvals
    p = [1 5 10 10 5 K];
    r = roots(p);
    plot(real(r), imag(r), 'x', 'DisplayName', ['K = ', num2str(K)]);
end

xline(0,'--k');
legend;

%% Problem 5
k = 5;
tau = 0.2;
s = tf('s');
G1 = k/(1 + tau*s);
figure; bode(G1); grid on


G1 = k/(1 + tau*s);
G2 = (s+2)/(s^2 + 3*s);
G3 = 20*(s+2)/(s^3 + 10*s^2 + 29*s + 20);

figure; bode(G2); grid on
figure; bode(G3); grid on

%% Problem 6

s = tf('s');

G = 20*(s+2)/((s+1)*(s^2 + 4*s + 100));

figure;
bode(G)
grid on