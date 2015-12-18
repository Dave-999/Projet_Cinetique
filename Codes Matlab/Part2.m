function [P] = Part2()

r_i = 8.36*10^-9;
AIBN=10^-5; %Considere comme constant
M0=0.5; %Valeur arbitraire plausibles
TrH0=10^-4;
f=0.7;
C_s=2.7; %Choix arbitraire d'un agent de tranfert

k_0=r_i/(2*f*AIBN);
k_t=10^7.7;
k_p=10^2.5;
k_s=C_s*k_p;

t=linspace(0,24*60*60,1000); %On analyse sur 24h

R=zeros(1,1000);
P=zeros(1,1000);
M=zeros(1,1000);
M(1)=M0;

alpha=zeros(1,1000);
TrH=zeros(1,1000);
TrH(1)=TrH0;

for i=1:1000
    R(i)=(r_i/k_t)^1/2;
    M(i)=M(1)*exp(-(r_i/k_t)^(1/2)*k_p*t(i));
    P(i)=r_i*t(i); 
    TrH(i)=TrH(1)*exp(-sqrt(r_i/k_t)*k_s*t(i));
    alpha(i)=k_p*M(i)/(k_p*M(i)+k_t*R(i));
end

%Les differents graphes obtenus
figure;
plot(t,ones(1,1000)*AIBN);
title('[AIBN]');
xlabel('t(s)');
ylabel('Concentration(1/L)');
axis([0 24*60*60 0 10^-4]);

figure;
plot(t,TrH);
title('[TrH]');
xlabel('t(s)');
ylabel('Concentration(1/L)');

figure;
plot(t,R);
title('[R]');
xlabel('t(s)');
ylabel('Concentration(1/L)');
axis([0 24*60*60 0 10^-15]);

figure;
plot(t,M);
title('[M]');
xlabel('t(s)');
ylabel('Concentration(1/L)');

figure;
plot(t,2./(1-alpha),t,(2+alpha)./(1-alpha));
title('DP');
legend('DPn','DPw');
xlabel('t(s)');

figure;
plot(t,P);
title('[PMMA]');
xlabel('t(s)');
ylabel('Concentration(1/L)');

end