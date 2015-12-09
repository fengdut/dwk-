
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



