#include<cmath>
#include<complex>

#include"tokamak.h"
#include"mconf.h"
#include"protos.h"
#include"AllocArray.h"
#include"simpintegral.h"
#include"vector.h"

using namespace std;
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

void Yps(const Grid *grid, complex<double> ** G_2D, double ** Chi_2D, double *** b_lambda_3D, double *** lambda_b_3D,double *** Theta_3D,int p,std::complex<double>** Yps_2D)
{
	using namespace std;
	complex<double> *tY;	
	Alloc1D(tY,grid->ntheta);
//	cout<<"Chi_2D ";
//	max_min_2D(grid->nx,grid->nL,Chi_2D);
//	cout<<"G_2D ";
//	max_min_2D(grid->nx,grid->nL,G_2D);
	
	for(int ix=0;ix<grid->nx;ix++)
		for(int iL=0;iL<grid->nL;iL++)
		{
			complex<double> texp=0;
			for(int it=0;it<grid->ntheta;it++)
			{
				complex<double> ti(0,-1.0);
				
				texp =exp(ti *Chi_2D[ix][iL]*(double)p*Theta_3D[ix][iL][it]);	
			//	texp=1.0;
				tY[it] =G_2D[ix][it] *b_lambda_3D[ix][iL][it] * (lambda_b_3D[ix][iL][it] +2*(1-lambda_b_3D[ix][iL][it]))*texp;
 
			}
			Yps_2D[ix][iL] = simpintegral(tY,grid->ntheta,grid->dtheta)*Chi_2D[ix][iL]/(2*M_PI);
		}	
	Free1D(tY);
}

void Yp_2( std::complex<double> ** const Yps_2D, Grid * const grid, double **Yp2)
{
	for(int ix=0;ix<grid->nx;ix++)
	for(int iL=0;iL<grid->nL;iL++)
	{
		Yp2[ix][iL] =abs(Yps_2D[ix][iL])*abs(Yps_2D[ix][iL]);	
	}
}

void omega_b(Grid * const grid, Tokamak * const tok,double ** const kappa, double ** const K, double * const q_1D,double *** omega_b_3D)
{
	for(int ix=0;ix<grid->nx;ix++)
	{
		double q=q_1D[ix];
		double x=grid->xarray[ix];
		for(int iL=0;iL<grid->nL;iL++)
		{
			double L =grid->Larray[iL];
			for(int iE=0;iE<grid->nE;iE++)
			{
				double E = grid->Earray[iE];
				omega_b_3D[ix][iL][iE] = M_PI *sqrt(kappa[ix][iL])/K[ix][iL]*sqrt(x*tok->eps*L*0.5)/q *sqrt(E);	
				//omega_b_3D[ix][iL][iE] = 1;	
			}
		}
	}

}

void tau_b(Grid *const grid, double *** const omega_b_3D, double *** tau_b_3D)
{
	for(int ix=0;ix<grid->nx;ix++)
        for(int iL=0;iL<grid->nL;iL++)
        for(int iE=0;iE<grid->nE;iE++)
	{
		//tau_b_3D[ix][iL][iE] = 2*M_PI/omega_b_3D[ix][iL][iE];
		tau_b_3D[ix][iL][iE] = 2*M_PI/omega_b_3D[ix][iL][iE];
	}

}

void omega_phi(Grid *const grid,double *const q_1D, double ***const omega_b_3D,double *** omega_phi_3D)
{
	for(int ix=0;ix<grid->nx;ix++)
	for(int iL=0;iL<grid->nL;iL++)
	for(int iE=0;iE<grid->nE;iE++)
	{
		omega_phi_3D[ix][iL][iE] = q_1D[ix] * omega_b_3D[ix][iL][iE];
//		omega_phi_3D[ix][iL][iE] = sqrt(grid->Earray[iE]);
	
	}	

}




