function [] = Part3_4()

r_i0 = 8.36*10^-9;
AIBN_0=10^-4;
f=0.7;
k_0=r_i0/(2*AIBN_0*f);
M0=0.5;
TrH0=10^-4;
C_s=2.7;

t=linspace(0,24*60*60,1000);

for i=1:1000
    
end
%M=zeros(1,1000);
%TrH=zeros(1,1000);
A=zeros(1000,2);
R=zeros(1,1000);
alpha=zeros(1,1000);
k_s=zeros(1,1000);
r_i=zeros(1,1000);
AIBN=zeros(1,1000);

k_t=zeros(1,1000);
X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
A=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];
p=polyfit(X,A,1);

k_p=zeros(1,1000);
V=[0.2 0.3 0.4 0.5 0.6 0.7 0.8];
W=[2.5 2.5 2.5 2.3 1.8 1.2 0.15];
q=polyfit(V,W,2);

[t,A]=ode45(@getA,t,[M0,TrH0]);

for i=1:1000
    M(i)=A(i,1);
    AIBN(i)=AIBN_0*exp(-k_0*t(i));
    r_i(i)=2*f*k_0*AIBN(i);
    k_t(i)=10.^(p(1)*(1-A(i,1)/M0)+p(2));
    if A(i,1) <= 0.4
        k_p(i)=10^2.5;
    elseif A(i,1) >= 0.8
        k_p(i)=1;
    else
        k_p(i)=10^(q(1)*(1-A(i,1)/M0)^2+q(2)*(1-A(i,1)/M0)+q(3));
    end
    R(i)=sqrt(r_i(i)/k_t(i));
    k_s(i)=C_s*k_p(i);
    alpha(i)=k_p(i)*A(i,1)/(k_p(i)*A(i,1)+k_t(i)*R(i)+k_s(i)*R(i)*A(i,2));
end

%Les différents graphes obtenus
figure;
plot(t,AIBN);
title('[AIBN]');
xlabel('t [s]');
ylabel('Concentration [mol/L]');

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