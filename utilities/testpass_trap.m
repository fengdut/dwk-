clear;
eps=0.4/1.65;
clf;

n=100;
L=linspace(0,1,n);

k=(1-L*(1-eps))./(2*eps*L);
plot(L,k);

hold all;

A=linspace(1,1,n);
plot(L,A);

ylim([0 2]);