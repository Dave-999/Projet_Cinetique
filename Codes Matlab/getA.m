function dAdt = getA(t,A)

r_i0 = 8.36*10^-9;
AIBN_0=10^-4;
AIBN=zeros(1,1000);
f=0.7;
k_0=r_i0/(2*AIBN_0*f);
M0=0.5;
TrH0=10^-4;
C_s=2.7;

r_i=2*f*k_0*AIBN_0*exp(-k_0*t);

k_t=zeros(1,1000);
X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
Y=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];
p=polyfit(X,Y,1);
k_t=10.^(p(1)*(1-A(1)/M0)+p(2));

k_p=zeros(1,1000);
V=[0.2 0.3 0.4 0.5 0.6 0.7 0.8];
W=[2.5 2.5 2.5 2.3 1.8 1.2 0.15];
q=polyfit(V,W,2)
k_p=10^(q(1)*(1-A(1)/M0)^2+q(2)*(1-A(1)/M0)+q(3));

k_s=C_s*k_p;

dAdt=zeros(2,1);
dAdt(1)=-k_p*A(1)*sqrt(r_i/k_t);
dAdt(2)=-k_s*A(2)*sqrt(r_i/k_t);

end