function [] = Part3_1()

r_i = 8.36*10^-9;
AIBN=10^-4; 
M0=0.5;
f=0.7;
C_s=2.7;

k_0=r_i/(2*f*AIBN);
k_p=10^2.5;
k_s=C_s*k_p;

t=linspace(0,24*60*60,1000);

M=zeros(1,1000);
R=zeros(1,1000);
alpha=zeros(1,1000);
k_t=zeros(1,1000);
X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
Y=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];
p=polyfit(X,Y,1);

[t,M]=ode45(@getM1,t,M0); 

for i=1:1000
    k_t(i)=10.^(p(1)*(1-M(i)/M0)+p(2));
    R(i)=sqrt(r_i/k_t(i));
    alpha(i)=k_p*M(i)/(k_p*M(i)+k_t(i)*R(i));
end

%Les différents graphes obtenus
figure;
plot(t,ones(1,1000)*AIBN);
title('[AIBN]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
axis([0 24*60*60 0 10^-3]);

figure;
plot(t,M);
title('[M]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');

figure;
plot(t,R);
title('[R]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
%axis([0 24*60*60]);

figure;
plot(t,2./(1-alpha),t,(2+alpha)./(1-alpha));
title('DP');
legend('DPn','DPw');
xlabel('t [s]');

%figure;
%plot(t,alpha);
%title('alpha');
%xlabel('t [s]');

end