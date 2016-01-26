clear;
clf(1);
figure(1);


data=load('dwk_omega_ED01.out');


romega=data(:,1);
iomega=data(:,2);
rdw=data(:,3);
idw=data(:,4);



n=1000;
omega_0=0,
omega_1=1.5;
omega=linspace(omega_0,omega_1,n)+0.002i;

dwk_kd=(1/3 +1/2*omega +omega.^2+omega.^3 .*log(1-1./omega))*20;



data=load('dwk_omega_ED001.out');




C=30;

 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 

plot(romega,rdw*C,'ko--');
hold all;

plot(omega,real(dwk_kd));





grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');


hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
plot(romega,idw*C,'ko--'); 
hold all;
plot(omega,imag(dwk_kd));


xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;
legend('$dwk++$','$analytical$');
title('$\omega_\star~ term$');

myprint('omega_dwk_omega_star_only_comp');


