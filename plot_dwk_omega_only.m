clear;
figure(1);

clf(1);
figure(1);




data=load('dwk_omega_only.out');
data1=load('dwk_omega_only.out1');

romega1=data1(:,1);
iomega1=data1(:,2);
rdw1=data1(:,3);
idw1=data1(:,4);


romega=data(:,1);
iomega=data(:,2);
rdw=data(:,3);
idw=data(:,4);



C=30;

 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.1 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 


hold all;
%plot(romega1,rdw1*C,'k+--');
plot(romega,rdw*C,'ro--');


%xlim([0 0.99]);

grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');


hax=axes('Position',[0.58 0.15 0.4 0.75],'FontSize',24); 
%xlim([0 0.99]);

hold all;
%plot(romega1,idw1*C,'k+--'); 

plot(romega,idw*C,'ro--'); 

%plot(romega,romega,'b--');
xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;


%myprint('omega_dwk');


