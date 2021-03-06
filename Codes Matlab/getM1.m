function dMdt = getM1(t,M)

r_i = 8.36*10^-9;
M0=0.5; 
k_p=10^2.5;

k_t=zeros(1,1000);
X=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1];
Y=[7.7 7.5 6.7 5.9 5.8 5.25 4.45 3.95 3.2 2];
p=polyfit(X,Y,1);
k_t=10.^(p(1)*(1-M/M0)+p(2));

dMdt = -k_p*M*sqrt(r_i/k_t);

end