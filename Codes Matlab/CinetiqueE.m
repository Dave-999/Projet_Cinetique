function [P] = CinetiqueE()
%Calcule l'évolution des concentrations en prenant en compte l'épuisement
%de l'AIBN

r_i1 = 8.36*10^-9;
AIBN0=10^-5;
f=0.7;
k_0=r_i1/(2*AIBN0*f); %k_0 considéré constant
M0=0.5;%Valeurs arbitraires plausibles
TrH0=10^-4;%Augmenter [TRH]t=0 diminue le DPn --> OK
t=linspace(0,24*60*60,1000);%On analyse sur 24h

AIBN=zeros(1,1000);
AIBN(1)=AIBN0;

for i=2:1000

    AIBN(i)=AIBN(1)*exp(-k_0*t(i)); %Epuisement de l'AIBN
 if AIBN(i)<=0
        break
    end
end
figure;
plot(t,AIBN);
%Fonction de kt
Xt=[0 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80];
kt=[10^7.7 10^7.5 10^6.7 10^5.9 10^5.8 10^5.25 10^4.45 10^3.95 10^3.2];
Pkt=polyfit(Xt,kt,1);

%Fonction de kp
Xp=[0.40 0.50 0.60 0.70 0.80];
kp=[10^2.5 10^2.3 10^1.8 10^1.2 10^0.15];
Pkp=polyfit(Xp,kp,2);

C_s=2.7; %Choix arbitraire d'un agent de tranfert

R=zeros(1,1000);
TrH=zeros(1,1000);
TrH(1)=TrH0;
P=zeros(1,1000);
M=zeros(1,1000);
M(1)=M0;

%R(0)=(r_i/10^7.7)^1/2;
alpha=zeros(1,1000);
for i=1:1000
    r_i=2*f*k_0*AIBN(i);
    k_t=polyval(Pkt,P(i));
    if P(i)<0.4
        k_p=10^2.5;
    end
    if P(i)>0.8
        k_p=0;
    end
    if P(i)>0.4 && P(i)<0.8
        k_p=polyval(Pkp,P(i));
    end
    R(i)=(r_i/k_t)^1/2;
    k_s=C_s*k_p;

    TrH(i)=TrH(1)*exp((2*k_s/k_0)*(R(i)-R(1))); % Cas général
    M(i)=M(1)*exp((2*k_p/k_0)*(R(i)-R(1)));
    P(i)=2*f*AIBN(1)*(1-exp(-k_0*t(i)))+k_s*R(1)*(exp((k_s/k_0)*R(1)*(exp(-k_0*t(i)/2)-1)-(k_0*t(i)/2))-1);
    alpha(i)=k_p*R(i)*M(i)/(k_p*R(i)*M(i)+k_s*R(i)*TrH(i)+k_t*R(i)^2);
end
figure
plot(t,2./(1-alpha))
title('DPn')
figure
plot(t,P)
title('[PMMA]')
figure
plot(t,M)
title('[M]')
end

