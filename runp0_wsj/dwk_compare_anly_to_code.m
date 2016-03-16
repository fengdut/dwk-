clf;
clear;

a=0.38;
R0=1.3;
p0=2.8083e3;
eps0=a/R0;

xi_s=0.01;
r_s=0.38;

E0=25.142e3;
B0=0.84;
mp=2;
Z=1;

rho_h=3.85698037e-02;
rho_h=rho_h/a;

delta_r=0.2;
c0=0.8;
c2=1.385;

n=100;

omegai=0.005i;
omega=linspace(0.1+omegai,0.95+omegai,n);

Ckd=-8*pi^2*a^2*R0*p0 *rho_h/eps0 *(eps0*xi_s)^2;
fs=exp (-(r_s/delta_r)^2)*(c2*(delta_r^2+r_s^2)+c0);
f0=c2*delta_r^2+c0;
df=fs-f0;
Ckd=Ckd*df;
omega_kd=omega.^3 .*log(1-1./omega) +omega.^2 +0.5*omega +1/3;
dwk_kd=omega_kd*Ckd;

Cks=4*pi^2*a^2*R0*p0*(eps0*xi_s)^2*delta_r^2;
dfs=exp(-r_s^2/delta_r^2) -1;
Cks=Cks*dfs;
omega_ks=-omega./(omega-1) +omega +omega.^2 .*log(1-1./omega);
dwk_ks =Cks *omega_ks;


dwk_k=dwk_kd+dwk_ks;
set(gcf,'Units','points','position',[100 500 1600 800],'Color',[1 1 1]);
hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 

plot(real(omega),real(dwk_kd),'-','LineWidth',2);
hold all;

plot(real(omega),real(dwk_ks),'g-','LineWidth',2);
grid on;

plot(real(omega),real(dwk_k),'r-','LineWidth',2);

xlabel('$real(\Omega)$');
ylabel('$real(\delta_W) ~(J)$');


%%%%read dwk++ data %%%%


Ct=0.47073675;


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_kd');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_rdwk_si,'bx');


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_ks');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_rdwk_si,'gx');


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_all');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_rdwk_si,'rx');

legend('$\delta W_{k,d}~(Analytical)~~$','$\delta W_{k,s}~(Analytical)~~$','$\delta W_{k}~(Analytical)~~$','$\delta W_{k,d} ~(dwk++)~$','$\delta W_{k,s} ~(dwk++)$','$\delta W_{k}~ ~(dwk++)$','Location','southwest');
hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 


grid on;

plot(real(omega),imag(dwk_kd),'-','LineWidth',2);
hold all;

plot(real(omega),imag(dwk_ks),'g-','LineWidth',2);
grid on;

plot(real(omega),imag(dwk_k),'r-','LineWidth',2);


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_kd');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_idwk_si,'bx');


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_ks');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_idwk_si,'gx');


[dwk_romega,dwk_iomega,dwk_rdwk,dwk_idwk]=read_omega_dwk('omega_dwk.out_all');
T0k=E0 *1.602e-19;
n0=p0/T0k;
n0=n0/Ct;
Cdwk=pi^2 *a^2 *R0 *n0*T0k
dwk_rdwk_si=dwk_rdwk*Cdwk;
dwk_idwk_si=dwk_idwk*Cdwk;
plot(dwk_romega,dwk_idwk_si,'rx');


legend('$\delta W_{k,d} ~(Analytical)~~$','$\delta W_{k,s}~(Analytical)~~$','$\delta W_{k}~(Analytical)~~$','$\delta W_{k,d} ~(dwk++)~$','$\delta W_{k,s} ~(dwk++)$','$\delta W_{k}~(dwk++)$','Location','southwest');

xlabel('$real(\Omega)$');
ylabel('$imag(\delta_W) ~(J)$');


%myprint('dwk_compare_anly_to_code');
