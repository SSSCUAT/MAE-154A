%% Homework 7
%% Adela Reyes 
clear, clc;

%% Problem 1

% Part 2 

T = 2.25;
c0 = -0.75;
w0 = (2*pi) / T;
t = -3:0.001:3;
k = -200:1:200;
u = c0 / T *ones(size(t)); 
w = k * w0;
j = 1i;

c_k = ((3 ./ (-j*k*w0)) .* exp(j*k*w0*0.25) .* (exp(-j*k*w0*1.25)- 1));
c_k(201) = c0;

u_FS = zeros(size(t));

for i = 1:length(w)
    u_FS = u_FS + (1/T)*(c_k(i)*exp(j*w(i)*t));
end

figure(1);
plot(t, real(u_FS));
grid on;
title('partial sum plot (u_FS)');
xlabel('time (s)');
ylabel('Output');
