clear;
figure(1);
clf(1);
figure(1);
B=0.84*1e4;

PB=3.98e6*(B/1e4)^2;

mh=2*1.67e-24;
xi0=1;

R0=1300;
a=38;



v0=190e3*2*pi*130;
E0=0.5*v0^2;
omega_0=190e3;

omega_c=9.58e3 *0.5 *B;


p0=0.01*PB;  %???
dr=0.3 *a;
rs=0.5*a;

Call=2*pi*R0*xi0^2*(rs*B/2/R0)^2;


C=4*2^1.5*B*mh*pi*pi*xi0*xi0/R0;
C1=-1*pi *R0*v0^3/omega_c;

C2= 1/(2^1.5 *pi*mh*B*E0);

%intr1=dr *exp(-rs^2/dr^2)*(2*dr^2+2*rs^2+0.5) -dr*(2*dr^2+0.5);
intr1=0.5*exp(-rs^2/dr^2) *(4*dr^2 +a^2+4*rs^2)/a^2 - 0.5*(4*dr^2+a^2)/a^2;
dwkd=C*C1*C2*intr1*p0;

Cs1=-pi*v0^2;
intr2=-0.5*dr^2*exp(-rs^2/dr^2) +0.5*dr^2;
dwks=C*Cs1*C2*intr2*p0;


n=100;
omega_a=0.1+0.01i;
omega_b=0.9+0.01i;
omega=linspace(omega_a,omega_b,n);

dwkdomega=omega.^3 .*log(1-1./omega) +omega.^2 +0.5.*omega +1/3;
dwksomega=-omega./(omega-1) +omega +omega.^2.*log(1-1./omega);

dwkdr=real(dwkd.*dwkdomega + dwks*dwksomega);

set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);

hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
 
 %rdwk=real(dwks*dwksomega)+real(dwkd.*dwkdomega);

rdwk=real(dwks.*dwkdomega);
plot(omega,real(dwks.*dwksomega)/1e7,'-','LineWidth',2);
grid on;

idwk=imag(dwks*dwksomega)+imag(dwkd.*dwkdomega);
idwk=idwk/Call;

hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
% plot(omega,dd_ds,'o--','LineWidth',2);
%plot(omega,imag(dwkd.*dwkdomega),'-','LineWidth',2);
plot(omega,imag(dwks*dwksomega)/1e7,'-','LineWidth',2);
grid on;