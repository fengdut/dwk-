clear;
 figure(2);
 clf;
 figure(2);
data=load('omega_dwk.out_star');


romega=data(:,1);
iomega=data(:,2);
rdwk=data(:,3);
idwk=data(:,4);

data=load('omega_dwk.out_omega');


romega1=data(:,1);
iomega1=data(:,2);
rdwk1=data(:,3);
idwk1=data(:,4);

set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);

hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 
plot(romega,rdwk,'ro--');
 
grid on;
 

xlabel('$real(\omega)$');
ylabel('$real(\delta W_k)$');




 hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
  plot(romega,rdwk1./rdwk,'ro--');


  grid on;
  xlabel('$real(\omega)$');
ylabel('$imag(\delta W_k)$');
