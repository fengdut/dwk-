
#include"tokamak.h"
#include"cmath"

void qprofile(const int nx,const double *xarray, double *q_1D)
{
	for(int i=0;i<nx;i++)
		q_1D[i] = 0.5 +2*xarray[i]*xarray[i];
	
}


double bf(const Tokamak *tok,const double theta, const double r)
{
	return 1-tok->eps*r*cos(theta);
}

void J_q(Grid *const grid,double *const q_1D, double *J_q)
{
	for(int ix=0;ix<grid->nx;ix++)
	{
		J_q[ix] = grid->xarray[ix]/q_1D[ix];
	}	
}



