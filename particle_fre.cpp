#include"tokamak.h"
#include"mconf.h"
#include"protos.h"
#include<cmath>
#include<complex>
#include"AllocArray.h"
#include"simpintegral.h"
void Chi(const Grid *grid,const Tokamak *tok, double sigma,double **Chi_2D,double **kappa_2D,double **K_2D)
{
	using namespace std;
	for(int ix=0;ix<grid->nx;ix++)
	{
		double teps =tok->eps *grid->xarray[ix];
		for(int iL=0;iL<grid->nL;iL++)
		{
			double k= (1-grid->Larray[iL] *(1-teps));
			k = k/(2*grid->Larray[iL]*teps);
			kappa_2D[ix][iL] = k;
			if(k>1) // for passing particles
			{
				double K=ellpk(1/k);
				K_2D[ix][iL]=K;
				Chi_2D[ix][iL] = sigma *M_PI *sqrt(k*teps *grid->Larray[iL]*0.5)/K;
			}
		}
	}
}

void Yps(const Grid *grid, double ** G_2D, double ** Chi_2D, double *** b_lambda_3D, double *** lambda_b_3D,double *** Theta_3D,int p,std::complex<double>** Yps_2D)
{
	using namespace std;
	complex<double> *tY;	
	Alloc1D(tY,grid->ntheta);
	for(int ix=1;ix<grid->nx;ix++)
		for(int iL=0;iL<grid->nL;iL++)
		{
			complex<double> texp;
			for(int it=0;it<grid->ntheta;it++)
			{
				complex<double> ti(0,-1.0);
				texp =exp(ti *Chi_2D[ix][iL]*(double)p*Theta_3D[ix][iL][it]);	
				tY[it] =G_2D[ix][it] *b_lambda_3D[ix][iL][it] * (lambda_b_3D[ix][iL][it] +2*(1-lambda_b_3D[ix][iL][it]))*texp;
 
			}
			Yps_2D[ix][iL] = simpintegral(tY,grid->ntheta,grid->dtheta)*Chi_2D[ix][iL]/(2*M_PI);
		}	
	Free1D(tY);
}





