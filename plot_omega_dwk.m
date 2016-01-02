clear;

ni=20;
nr=20;

data=load('dwk.out');
romega=data(:,1);
iomega=data(:,2);

rdwk=data(:,3);
idwk=data(:,4);

romega=reshape(romega,[ni,nr]);
iomega=reshape(iomega,[ni,nr]);
rdwk=reshape(rdwk,[ni,nr]);

idwk=reshape(idwk,[ni,nr]);


[X,Y] = meshgrid(romega(:,1),iomega(1,:));
subplot(1,2,1);
contourf(X',Y',rdwk,50,'LineStyle','none');
colorbar;
colormap(jet);
subplot(1,2,2);
contourf(X',Y',idwk,50,'LineStyle','none');
colorbar;
colormap(jet);
