function [P] = Part2()

r_i = 8.36*10^-9;
AIBN=10^-5; 
M0=0.5; 
TrH0=10^-4;
f=0.7;
C_s=2.7;

k_0=r_i/(2*f*AIBN);
k_t=10^7.7;
k_p=10^2.5;
k_s=C_s*k_p;

R=sqrt(r_i/k_t);

t=linspace(0,24*60*60,1000); %On analyse sur 24h

M=ones(1,1000);
TrH=ones(1,1000);
alpha=ones(1,1000);

for i=1:1000
    M(i)=M0*exp(-R*k_p*t(i)); 
    TrH(i)=TrH0*exp(-R*k_s*t(i));
    alpha(i)=k_p*M(i)/(k_p*M(i)+k_t*R);
end

%Les differents graphes obtenus
figure;
plot(t,ones(1,1000)*AIBN);
title('[AIBN]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
axis([0 24*60*60 0 10^-4]);

figure;
plot(t,M);
title('[M]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');

figure;
plot(t,ones(1,1000)*R);
title('[R]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
axis([0 24*60*60 10^-8 10^-7]);

figure;
plot(t,TrH);
title('[TrH]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');

figure;
plot(t,2./(1-alpha),t,(2+alpha)./(1-alpha));
title('DP');
legend('DPn','DPw');
xlabel('t [s]');

end