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
c2=1.3842;

% c0=1;
% c2=0;

epss=r_s*eps0;

dr=0.2;
omega_r=0.8472;
%omega_r=0.8901;

omega_A=2.935;
A=dr*dr *(exp(-(r_s/dr)^2)- 1);%
%A=0;


B= 2*rho_h /eps0 *omega_r * (exp(-(r_s/dr)^2)*(c2*(dr^2+r_s^2)+c0) -(c2*(dr^2)+c0) );

beta_crit=epss^2 /(pi *omega_A*omega_r *(A-B)*eps0*eps0)