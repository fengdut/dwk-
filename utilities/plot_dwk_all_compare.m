clear;
figure(1);
clf(1);
figure(1);


data=load('dwk_omega_dwk.out');

romega=data(:,1);
iomega=data(:,2);
rdw=data(:,3);
idw=data(:,4);

n=1000;
omega_0=0,
omega_1=1.1;
omega=linspace(omega_0,omega_1,n)+0.002i;

dwk_omega=-0.002*real(omega).*((omega-2)./(omega-1)+omega.*log(1-1./omega));
dwk_omegastar=0.0205*(1/3 +1/2*omega +omega.^2+omega.^3 .*log(1-1./omega));
dwk_all=dwk_omega+dwk_omegastar;

set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
hax=axes('Position',[0.1 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 

plot(real(omega),real(dwk_all),'k.','LineWidth',4);
hold all;
plot(romega,rdw,'ro--','MarkerSize',8);
ylim([0-0.12 0.05]);

xlim([0 1.1]);

grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

legend('$analytical~~$','$dwk++$','Location','northwest');
text(0.4,0.06,'$dwk++:~\Lambda_0=0.01,\Delta \Lambda=0.01,E_d=0.01, E_c=0.01, \Delta r=0.001$','FontSize',24);

hax=axes('Position',[0.58 0.15 0.4 0.75],'FontSize',24); 
hold all;

plot(real(omega),imag(dwk_all),'k.','LineWidth',4);
plot(romega,idw,'ro--','MarkerSize',8); 


ylim([0-0.12 0.05]);
xlim([0 1.1]);
xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;

myprint('omega_dwk++_analytical_l0.01_ld0.01_ed0.01_ec0.01_dr0.001_all');


