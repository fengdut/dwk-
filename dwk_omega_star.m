clear;
clf(2);
figure(2);

n=1000;
omega_0=0,
omega_1=1.5;
omega=linspace(omega_0,omega_1,n);

dwk_kd=1/3 +1/2*omega +omega.^2+omega.^3 .*log(1-1./omega);

 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
plot(omega,real(dwk_kd));


grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
plot(omega,imag(dwk_kd));

xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;


myprint('wangsj_dwk_omega_star');