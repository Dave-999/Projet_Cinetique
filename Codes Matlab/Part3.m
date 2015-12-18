function [] = Part3()

close all;
r_i = 8.36*10^-9;
AIBN0=10^-5;
M0=0.5; %Valeurs arbitraires plausibles
TrH0=10^-4; %Augmenter [TRH]t=0 diminue le DPn --> OK
f=0.7;
t=linspace(0,24*60*60,1000); %On analyse sur 24h

AIBN=AIBN0;


%Fonction de kt
%Xt=[0 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80];
Xt = [0.10 1]
%kt=[10^7.7 10^7.5 10^6.7 10^5.9 10^5.8 10^5.25 10^4.45 10^3.95 10^3.2];
kt= [10^7.5 10^2]
Pkt=polyfit(Xt,kt,1)
%figure
%plot(linspace(0,1,1000),polyval(Pkt,linspace(0,1,1000)))

%Fonction de kp
Xp=[0.40 0.50 0.60 0.70 0.80];
kp=[10^2.5 10^2.3 10^1.8 10^1.2 10^0.15];
Pkp=polyfit(Xp,kp,2)

C_s=2.7; %Choix arbitraire d'un agent de tranfert

R=zeros(1,1000);
TrH=zeros(1,1000);
TrH(1)=TrH0;
P=zeros(1,1000);
M=zeros(1,1000);
M(1)=M0;
k_p=10^2.5;
%R(0)=(r_i/10^7.7)^1/2;
alpha=zeros(1,1000);
for i=2:1000

    %3 EQUATIONS A 3 INCONNUES QUE JE N4ARRIVE PAS A RESOUDRE!!!!
    %XM=(1-M(i)/M(1));
    %int1=2*sqrt(Pkt(2)+Pkt(1)*XM)/(Pkt(1))-2*sqrt(Pkt(2))/(Pkt(1)) %Intégrale de kt^-1/2
    %M(i)=M(1)*exp(-(r_i^(1/2)*k_p*int1));
    
    k_t=polyval(Pkt,(1-M(i)/M(1)));
    R(i)=(r_i/k_t)^1/2;
    %k_s=C_s*k_p; %A décommenter si agent de transfert
    %TrH(i)=TrH(1)*exp(-sqrt(r_i/k_t)*k_s*t(i)); %Si r_i = cste
    P(i)=r_i*t(i);%-(k_s/k_t)*exp(-sqrt(r_i/k_t)*k_t*t(i)) + k_s/k_t;
    alpha(i)=k_p*R(i)*M(i)/(k_p*R(i)*M(i)+k_t*R(i)^2);%+k_s*R(i)*TrH(i));
end


%Les différents graphes obtenus
figure
plot(t,ones(1,1000)*AIBN)
title('[AIBN]')
xlabel('t(s)')
ylabel('Concentration(1/L)')
axis([0 24*60*60 0 10^-4]);

figure
plot(t,R)
title('[R]')
xlabel('t(s)')
ylabel('Concentration(1/L)')
axis([0 24*60*60 0 10^-15]);

figure
plot(t,M)
title('[M]')
xlabel('t(s)')
ylabel('Concentration(1/L)')

figure
plot(t,2./(1-alpha),t,(2+alpha)./(1-alpha))
title('DP')
legend('DPn','DPw')
xlabel('t(s)')

figure
plot(t,P)
title('[PMMA]')
xlabel('t(s)')
ylabel('Concentration(1/L)')

end