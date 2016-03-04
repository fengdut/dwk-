clear;

R=13.0;
a=0.38;
xi=0.01;

Th=21.90995e3;



dwkd=0.001187519886;
dwkd=0.002776703252;

dwks=-0.000030387092;



Cn= 0.540278880838;
CP= 0.278328331929;


n=8e17;
n=n*Cn/CP;
C=pi*pi*xi*xi*R *n*Th*1.602176e-19

dwkd =dwkd*C
dwks =dwks*C


figure(2);
clf(2);
figure(2);
data=load('omega_dwk.out_dwks');
%data=load('test.out');

romega=data(:,1);
iomega=data(:,2);

rdwk=data(:,3);
idwk=data(:,4);

 set(gcf,'Units','points','position',[10 80 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
 plot(romega,rdwk*C,'ro--');
 
 grid on;
 
%  ylim([0-0.12 0.1]);
% xlim([0 1.1]);
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');

%text(0.4,0.11,'$dwk++:~\Lambda_0=0.28,\Delta \Lambda=0.1,\Delta E=0.3, E_c=0.25, \Delta r=0.001$','FontSize',24);



 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(romega,idwk*C,'ro--');

% ylim([0-0.12 0.15]);
% xlim([0 1.1]);
  grid on;
  xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');
  
 % myprint('dwk_omega_dwk');


