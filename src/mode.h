
#ifndef MODE_H
#define MODE_H
#include"tokamak.h"
#include<complex>

typedef struct MODE
{
	int n;
	int m;
	int pa,pb;
	double delta_r;
	
	double omega_0,omega_1;
	double omega_i;
	int omega_n;
	std::complex<double> *omega_array;
	
	double omega_err;
	int max_iter;
	int max_iterg;
	double dw_f;
}Mode;

//void G_R_theta(Grid * const grid, Tokamak * const tok,Mode * const pmode,double **G_2D);
void G_R_theta(Grid * const grid, Tokamak * const tok, Slowing *const slow,Mode * const pmode,double **G_2D);
void G_R_theta(Grid * const grid, Tokamak * const tok, Slowing *const slow,Mode * const pmode, double *const q_1D,std::complex<double> **G_2D);

#endif
