filename0='omega_dwk.out_p1';
filename1='omega_dwk.out_p01';
    
[aomegar,aomegai,adwkr,adwki]=read_omega_dwk('omega_dwk.out_p1');
[bomegar,bomegai,bdwkr,bdwki]=read_omega_dwk('omega_dwk.out_p01');
[comegar,comegai,cdwkr,cdwki]=read_omega_dwk('omega_dwk.out_p0');
[domegar,domegai,ddwkr,ddwki]=read_omega_dwk('omega_dwk.out_p-1');



 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
 plot(aomegar,adwkr,'ro--');
 
 hold all;

plot(bomegar,bdwkr,'bo--');
plot(comegar,cdwkr,'ko--');
plot(domegar,ddwkr,'go--');



 
 
 grid on;
 
%  ylim([0-0.12 0.1]);
% xlim([0 1.1]);
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

%text(0.4,0.11,'$dwk++:~\Lambda_0=0.28,\Delta \Lambda=0.1,\Delta E=0.3, E_c=0.25, \Delta r=0.001$','FontSize',24);



 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(aomegar,adwki,'ro--');
hold all;
    plot(bomegar,bdwki,'bo--');
        plot(comegar,cdwki,'ko--');

        plot(domegar,ddwki,'go--');

% ylim([0-0.12 0.15]);
% xlim([0 1.1]);
  grid on;
  xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');
  
 % myprint('dwk_omega_dwk');
 
