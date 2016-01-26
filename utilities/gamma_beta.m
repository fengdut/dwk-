clear;
clf;

gamma=[0.001            0.002            0.004           0.008          0.012               0.016           0.02            0.025                   0.03];
omega=[0.880530166626   0.880056095123   0.879152965546  0.877464962006 0.875923824310      0.874520969391  0.873249721527  0.871835422516        0.870606136322];
beta_h=[0.031930639361  0.032359585558   0.033210189630  0.034970529593 0.036814950984      0.038749959036  0.040782435211  0.043471629399        0.046340875013];


 set(gcf,'Units','points','position',[100 500 1200 600],'Color',[1 1 1]);
% 
% % 
 hax=axes('Position',[0.08 0.15 0.4 0.75],'FontSize',24,'FontName','Latex'); 

plot(beta_h,gamma,'ro--','LineWidth',2);
xlabel('$\beta_h$');
ylabel('$\gamma$');
ylim([0 0.035]);
grid on;
hax=axes('Position',[0.56 0.15 0.4 0.75],'FontSize',24); 
plot(beta_h,omega,'ro--','LineWidth',2);
xlabel('$\beta_h$');
ylabel('$\omega$');
grid on;

title('P=0');

myprint('test_gamma_beta');