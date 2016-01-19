clear;
clf;

data=load('dwk_omega_dwk.out');
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
 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(romega,idwk,'ro--');

  grid on;