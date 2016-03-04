clear;
clf;

ncid=netcdf.open('Yps.nc','NC_NOWRITE');
id=netcdf.inqVarID(ncid,'r');
[varname,xtype,varDminIDs,varAtts]=netcdf.inqVar(ncid,id);
r=netcdf.getVar(ncid,id);
nr=size(r,1);

id=netcdf.inqVarID(ncid,'L');
[varname,xtype,varDminIDs,varAtts]=netcdf.inqVar(ncid,id);
L=netcdf.getVar(ncid,id);
nL=size(L,1);

id=netcdf.inqVarID(ncid,'E');
[varname,xtype,varDminIDs,varAtts]=netcdf.inqVar(ncid,id);
E=netcdf.getVar(ncid,id);
nE=size(E,1);

id=netcdf.inqVarID(ncid,'rYps_-1');
[varname,xtype,varDminIDs,varAtts]=netcdf.inqVar(ncid,id);
rYps=netcdf.getVar(ncid,id);

id=netcdf.inqVarID(ncid,'iYps_-1');
[varname,xtype,varDminIDs,varAtts]=netcdf.inqVar(ncid,id);
iYps=netcdf.getVar(ncid,id);


netcdf.close(ncid);

iE=100;
iL=100;

rYps1(1:nr)=rYps(iE,iL,1:nr);



 set(gcf,'Units','points','position',[100 500 800 600],'Color',[1 1 1]);
 hax=axes('Position',[0.15 0.15 0.75 0.75],'FontSize',24,'FontName','Latex'); 

 
plot(r,rYps1(:),'r.--');
titstr=sprintf('$p=-1,E=%f,\\Lambda=%f$',E(iE),L(iL));
title(titstr);


r_s=0.349773148705;
dr=0.01;

nrs=10;
rsarray=linspace(r_s,r_s,nrs);
y=linspace(-0.02, 0.1,nrs);
hold all;
plot(rsarray,y,'k-','LineWidth',2);

xlim([0 1]);
grid on;
xlabel('$r$');
ylabel('$real(Yps)$');

%myprint('Yps_-1_E05');


