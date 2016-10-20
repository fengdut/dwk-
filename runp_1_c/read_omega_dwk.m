function [omegar,omegai,dwkr,dwki]=read_omega_dwk(filename)



data0=load(filename);

omegar=data0(:,1);
omegai=data0(:,2);

dwkr=data0(:,3);
dwki=data0(:,4);

end