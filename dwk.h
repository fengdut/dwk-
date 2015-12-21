#ifndef DWK_H
#define DWK_H

#include<complex>

using namespace std;

complex<double> dwk(Grid *const grid, complex<double> const omega, int n,int p,
        double *** const omega_phi_3D, double ***const omega_b_3D, double *** tau_b_3D,
        double ** const Yp2_2D, double * const J_q_1D, double *** const F_E_3D, double *** const omega_star,complex<double> *** dwk_3D );
void help();


#endif
