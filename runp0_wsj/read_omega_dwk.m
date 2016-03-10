function [romega,iomega,rdwk,idwk]=read_omega_dwk(filename)

if(nargin<1)
    filename='omega_dwk.out';
end


data=load(filename);

romega=data(:,1);
iomega=data(:,2);

rdwk=data(:,3);
idwk=data(:,4);
end
