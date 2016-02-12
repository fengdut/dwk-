clear;
 figure(2);
 clf;
 figure(2);
data=load('omega_dwk.out');
%data=load('test.out');

romega=data(:,1);
iomega=data(:,2);

rdwk=data(:,3);
idwk=data(:,4);

 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
 plot(romega,rdwk,'ro--');
 
 grid on;
 
%  ylim([0-0.12 0.1]);
% xlim([0 1.1]);
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

%text(0.4,0.11,'$dwk++:~\Lambda_0=0.28,\Delta \Lambda=0.1,\Delta E=0.3, E_c=0.25, \Delta r=0.001$','FontSize',24);



 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(romega,idwk,'ro--');

% ylim([0-0.12 0.15]);
% xlim([0 1.1]);
  grid on;
  xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');
  
 % myprint('dwk_omega_dwk');