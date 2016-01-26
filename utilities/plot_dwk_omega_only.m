clear;
figure(1);

clf(1);
figure(1);



data=load('dwk_omega_only.out');


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

plot(romega,rdw*C,'ro--');


xlim([0 0.99]);

grid on;
xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');


hax=axes('Position',[0.58 0.15 0.4 0.75],'FontSize',24); 
xlim([0 0.99]);

hold all;


plot(romega,idw*C,'ro--'); 


xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');

grid on;


myprint('omega_dwk_only_dwk++');


