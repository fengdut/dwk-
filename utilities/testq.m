clear;
clf;
r=[0    0.15 0.3     0.45 0.8 1];
%q=[0.95  0.95   0.98  1   2 2.9];
%q=[0.92 +2.5*r];

%plot(r,q,'o--');

%cq=[0.87 1.6 -8.185 15.82 -7.022];
cq=[0.92 0 2.5 0 0];

ra=linspace(0,1,100);
qa=cq(1) +cq(2).*ra +cq(3).*ra.^2 +cq(4).*ra.^3 +cq(5).*ra.^4;
hold all;
plot(ra,qa);
grid on;




