clear;
figure(2);
%clf(2);
figure(2);

n=1000;
omega_0=0,
omega_1=0.99;
omega=linspace(omega_0,omega_1,n);

%dwk_kd0=-(-omega./(omega-1)+omega+omega.^2.*log(1-1./omega));

dwk_kd0=-((omega.^2-2*omega)./(omega-1)+omega.^2.*log(1-1./omega));

%dwk_kd=omega.^2./(omega-1);



 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
%plot(omega,real(dwk_kd));
hold all;
plot(omega,real(dwk_kd0));


grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
%plot(omega,imag(dwk_kd));
plot(omega,imag(dwk_kd0));

xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;


%myprint('wangsj_dwk_omega_only');