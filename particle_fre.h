#ifndef PARTICLE_FRE_H
#define PARTICLE_FRE_H
#include<complex>

void Chi(const Grid *grid,const Tokamak *tok, double sigma,double **Chi_2D,double **kappa_2D,double **K_2D);
void Yps(const Grid *grid, double ** G_2D, double ** Chi_2D, double *** b_lambda_3D, double *** lambda_b_3D,double *** Theta_3D,int p,std::complex<double>** Yps_2D);

#endif


