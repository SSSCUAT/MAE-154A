% MAE 171A
% Homework 1
% Adela Reyes 

%% Problem 1

k1=[100, 500, 1000, 5000];
time=linspace(0,40);
vt=zeros(1,40);

figure(1)
title("v(t) with different k values")
hold on
grid on

for j=1:1:length(k1)
    k=k1(j);
    for i=1:1:length(time)
        t=time(i);
        v=(k/(k+50))*(1-exp(((k/1000)+(1/20))*-t));
        vt(i)=v;
    end
    plot(time, vt);
end
hold off
