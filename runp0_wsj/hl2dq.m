clear;
clf;

n=50;

psi=linspace(0,1,n);

ac=[0.905753017013633, 0.223336609969303,2.983896706575029,-10.371299037980963, 22.771789907875711,-22.537492570337179 ,9.371000135373542];
q=ac(1) +ac(2)*psi +ac(3) *psi.^2 +ac(4) *psi.^3 +ac(5) *psi.^4  +ac(6)*psi.^5 +ac(7)*psi.^6;


plot(psi,q,'k--');

r=psi.^0.5;
plot(r,q);

rc=[0.009053654243265   0.012296316022935  -0.185406103486829   1.085989546198082  -3.061547275476543   4.601983792960791  -3.526051652813747   1.097082614662107]*100;
rc=[0.008053654243265   0.01296316022935  -0.185406103486829   1.085989546198082  -3.061547275476543   4.601983792960791  -3.526051652813747   1.097082614662107]*100;


r=linspace(0,1,n);
qr=rc(1) + rc(2) *r +rc(3) *r.^2 +rc(4) *r.^3 +rc(5) *r.^4 +rc( 6)*r.^5 +rc(7) *r.^6 +rc(8)*r.^7;

hold on;

plot(r,qr,'r-');
grid on;