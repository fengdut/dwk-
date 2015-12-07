
#include"tokamak.h"

void qprofile(int nx,const double *xarray, double *q_1D)
{
	for(int i=0;i<nx;i++)
		q_1D[i] = 0.5 +2*xarray[i]*xarray[i];
	
}





