function [P] = Cinetique()

r_i = 8.36*10^-9;
AIBN0=10^-5;
M0=0.5;%Valeurs arbitraires plausibles
TrH0=10^-4;%Augmenter [TRH]t=0 diminue le DPn --> OK
f=0.7;
t=linspace(0,24*60*60,1000);%On analyse sur 24h

AIBN=zeros(1,1000)
AIBN(1)=AIBN0;

for i=2:1000
AIBN(i)=AIBN(1);%Attention, on considère r_i constante ici -> AIBN inépuisable 

    %    AIBN(i)=AIBN0-(r_i/(2*f))*t(i); 
 %   if AIBN(i)<=0
  %      break
   % end
end

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
    k_0=r_i/(2*f*AIBN(1));
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
    R(i)=(r_i/k_t)^1/2
    k_s=C_s*k_p;
    TrH(i)=TrH(1)*exp(-sqrt(r_i/k_t)*k_s*t(i)); %Si r_i = cste
    M(i)=M(1)*exp(-sqrt(r_i/k_t)*k_p*t(i));
    P(i)=r_i*t(i)-(k_s/k_t)*exp(-sqrt(r_i/k_t)*k_t*t(i))+k_s/k_t;
    %TrH(i)=TrH0-exp((2*k_s/k_0)*(R(i)-R(1))); % Cas général
    %M(i)=M0-exp((2*k_p/k_0)*(R(i)-R(1)));
    %P(i)=2*f*AIBN(1)*(1-exp(-k_0*t(i)))+k_s*R(1)*(exp((k_s/k_0)*R(1)*(exp(-k_0*t(i)/2)-1)-(k_0*t(i)/2))-1);
    alpha(i)=k_p*R(i)*M(i)/(k_p*R(i)*M(i)+k_s*R(i)*TrH(i)+k_t*R(i)^2)
end
figure
plot(t,2./(1-alpha))
title('DPn')
figure
plot(t,P)
title('[PMMA]')
end
