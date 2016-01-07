clear;
figure(1);

clf(1);
figure(1);
data=load('dwk_omega_only_dlde.out');

romega=data(:,1);
iomega=data(:,2);
rdw=data(:,3);
idw=data(:,4);



data=load('dwk_omega_only.out');

romega1=data(:,1);
iomega1=data(:,2);
rdw1=data(:,3);
idw1=data(:,4);



C=30;

 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.07 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 

plot(romega,rdw*C,'ko--');
hold all;
plot(romega1,rdw1*C,'r+--');


grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');


hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
plot(romega,idw*C,'ko--'); 
hold all;
plot(romega1,idw1*C,'r+--'); 


%plot(romega,romega,'b--');
xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;


%myprint('omega_dwk');


