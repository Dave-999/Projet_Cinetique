function [] = Part3bis()

r_i = 8.36*10^-9;
AIBN=10^-5; 
M0=0.5; 
TrH0=10^-4;
f=0.7;
C_s=2.7;

k_0=r_i/(2*f*AIBN);
k_p=10^2.5;
k_s=C_s*k_p;

t=linspace(0,24*60*60,1000);

M=zeros(1,1000);
k_t=zeros(1,1000);

X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
Y=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];

p=polyfit(X,Y,1);
k_t=10.^(p(1)*(1-M/M0)-p(2));

[t,M]=ode45(@getM,t,M0); 

figure;
plot(t,M);
title('[M]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
axis([0 24*60*60 0 1]);

end