clear;
 figure(1);
 clf;
 figure(1);
data=load('omega_dwk_dr03_p-1.out');

romega=data(:,1);
iomega=data(:,2);
rdwk=data(:,3);
idwk=data(:,4);

data=load('omega_dwk_dr03_p0.out');
romega1=data(:,1);
iomega1=data(:,2);
rdwk1=data(:,3);
idwk1=data(:,4);


data=load('omega_dwk_dr03_p0-1.out');
romega2=data(:,1);
iomega2=data(:,2);
rdwk2=data(:,3);
idwk2=data(:,4);


 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
 plot(romega,rdwk,'bo--');
 hold all;
 plot(romega1,rdwk1,'ko--');
 grid on;
 
 
  plot(romega2,rdwk2,'ro--');
  
%  ylim([0-0.12 0.1]);
% xlim([0 1.1]);
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

%text(0.4,0.11,'$dwk++:~\Lambda_0=0.28,\Delta \Lambda=0.1,\Delta E=0.3, E_c=0.25, \Delta r=0.001$','FontSize',24);



 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(romega,idwk,'bo--');
hold all;
 plot(romega1,idwk1,'ko--');
  plot(romega2,idwk2,'ro--');
% ylim([0-0.12 0.15]);
% xlim([0 1.1]);
  grid on;
  xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');
  
 myprint('dwk_omega_dwk_p0_p1_pt');