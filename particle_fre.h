#ifndef PARTICLE_FRE_H
#define PARTICLE_FRE_H
#include<complex>

void Chi(const Grid *grid,const Tokamak *tok, double sigma,double **Chi_2D,double **kappa_2D,double **K_2D);
void Yps(const Grid *grid, std::complex<double> ** G_2D, double ** Chi_2D, double *** b_lambda_3D, double *** lambda_b_3D,double *** Theta_3D,int p,std::complex<double>** Yps_2D);
void Yp_2( std::complex<double> ** const Yps_2D, Grid * const grid, double **Yp2);

void omega_b(Grid * const grid, Tokamak * const tok,double ** const kappa, double ** const K, double * const q_1D,double *** omega_b_3D);
void omega_phi(Grid *const grid,double *const q_1D, double ***const omega_b_3D,double *** omega_phi_3D);
void tau_b(Grid *const grid, double *** const omega_b_3D, double *** tau_b_3D);
#endif


