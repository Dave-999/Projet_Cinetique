function [] = Part3_3()

r_i = 8.36*10^-9;
AIBN=10^-4; 
M0=0.5;
TrH0=10^-4;
f=0.7;
C_s=0.63;

k_0=r_i/(2*f*AIBN);

t=linspace(0,24*60*60,1000);

%M=zeros(1,1000);
%TrH=zeros(1,1000);
Z=zeros(1000,2);
R=zeros(1,1000);
alpha=zeros(1,1000);
k_s=zeros(1,1000);

k_t=zeros(1,1000);
X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
Y=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];
p=polyfit(X,Y,1);

k_p=zeros(1,1000);
V=[0.2 0.3 0.4 0.5 0.6 0.7 0.8];
W=[2.5 2.5 2.5 2.3 1.8 1.2 0.15];
q=polyfit(V,W,2);

[t,Z]=ode45(@getZ,t,[M0,TrH0]);

for i=1:1000
    M(i)=Z(i,1);
    k_t(i)=10.^(p(1)*(1-Z(i,1)/M0)+p(2));
    if Z(i,1) <= 0.4
        k_p(i)=10^2.5;
    elseif Z(i,1) >= 0.8
        k_p(i)=1;
    else
        k_p(i)=10^(q(1)*(1-Z(i,1)/M0)^2+q(2)*(1-Z(i,1)/M0)+q(3));
    end
    R(i)=sqrt(r_i/k_t(i));
    k_s(i)=C_s*k_p(i);
    alpha(i)=k_p(i)*Z(i,1)/(k_p(i)*Z(i,1)+k_t(i)*R(i)+k_s(i)*R(i)*Z(i,2));
end

%Les différents graphes obtenus
figure;
plot(t,ones(1,1000)*AIBN);
title('[AIBN]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');
axis([0 24*60*60 0 10^-3]);

figure;
plot(t,1-M/M0);
title('X_M');
xlabel('t [s]');
ylabel('Fraction convertie en monomères');

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

% figure;
% plot(t,alpha);
% title('alpha');
% xlabel('t [s]');

end