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

dwk_kd0=-0.002*((omega.^2-2*omega)./(omega-1)+omega.^2.*log(1-1./omega));

%dwk_kd0=-0.006*(omega.*(omega.*log(1-1./omega)+1));  %a

%dwk_kd0=0.002*(omega./(omega-1));                  %b

%dwk_kd0=0.006*(omega.*(omega.*log(1-1./omega)+1))/1.5; %c

set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.1 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 


hold all;
%plot(romega1,rdw1*C,'k+--');
plot(romega,rdw,'ro--');
plot(omega,real(dwk_kd0));
ylim([0-0.25 0.25]);

%xlim([0 0.99]);

grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');


hax=axes('Position',[0.58 0.15 0.4 0.75],'FontSize',24); 
%xlim([0 0.99]);
hold all;
%plot(romega1,idw1*C,'k+--'); 

plot(romega,idw,'ro--'); 
plot(omega,imag(dwk_kd0));

ylim([0-0.25 0.25]);

%plot(romega,romega,'b--');
xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;
title('$omega ~term$');
%legend('$dwk++$','$analytical$');

%myprint('omega_dwk_only_dwk++_analytical_l0.01_ld0.05_ed0.05_ec0.01_dr0.001');


