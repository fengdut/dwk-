#include<cmath>
#include<complex>

#include"tokamak.h"
#include"mconf.h"
#include"protos.h"
#include"AllocArray.h"
#include"simpintegral.h"
#include"vector.h"
#include"omp.h"
#include<ctime>
#include"stdlib.h"
#include"outlog.h"

using namespace std;
void Chi(const Grid *grid,const Tokamak *tok, double sigma,double **Chi_2D,double **kappa_2D,double **K_2D)
{
	using namespace std;
#pragma omp parallel for
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
					K_2D[ix][iL]=K;  //passing particles
				Chi_2D[ix][iL] = sigma *M_PI *sqrt(k*teps *grid->Larray[iL]*0.5)/K;
			}
			else
			{
				K_2D[ix][iL]=ellpk(k);	//trapped particles
				cout<<"trapped"<<endl;
				cerr<<"trapped particle is not finish"<<endl;
				exit(-1);
			}	
		}
	}
}

void Yps(const Grid *grid, complex<double> *** G_3D, double ** Chi_2D, double *** b_lambda_3D, double *** lambda_b_3D,double *** Theta_3D,int p,std::complex<double>*** Yps_3D)
{
	using namespace std;
   	clock_t c_start=clock();
	
	#pragma omp parallel 
	{	
	complex<double>*tY;	
	complex<double>*tYb;	
	Alloc1D(tY,grid->ntheta);
	Alloc1D(tYb,grid->ntheta);
	int id=omp_get_thread_num();
	int np=omp_get_num_threads();

	int ni=grid->nx/np;
	int nm=0;
	if(id+1==np)
		nm=grid->nx%np;	

/*	 cout<<"G_3D:\t";
        max_min_3D(grid->nx,grid->nE,grid->ntheta,G_3D);

	  cout<<"b_lambda:\t";
        max_min_3D(grid->nx,grid->nL,grid->ntheta,b_lambda_3D);

	   cout<<"lambda_b:\t";
        max_min_3D(grid->nx,grid->nL,grid->ntheta,lambda_b_3D);
*/

	for(int ix=id*ni;ix<(id+1)*ni+nm;ix=ix+1)
	{
		complex<double> texp=0;
		complex<double> ti(0,-1.0);
		for(int iL=0;iL<grid->nL;iL++)
		{
			double tchi=Chi_2D[ix][iL];
			double tchi2pi=tchi/(2.0*M_PI);		
			
			for(int it=0;it<grid->ntheta;it++)
			{
				texp 	=exp(ti*tchi*(double)p*Theta_3D[ix][iL][it]);
                                tY[it] 	=b_lambda_3D[ix][iL][it] * (2.0-lambda_b_3D[ix][iL][it])*texp;		
			}
			for(int iE=0;iE<grid->nE;iE++)
			{
				for(int it=0;it<grid->ntheta;it++)
				{
					tYb[it] =G_3D[ix][iE][it] *(tY[it]);
				}
				Yps_3D[ix][iL][iE] = simpintegral_o_p(tYb,grid->ntheta,grid->dtheta)*tchi2pi;
			}
		}	
	}
	Free1D(tY);
	Free1D(tYb);
	}
//      cout<<"time:\t"<<float(clock() - c_start)/CLOCKS_PER_SEC<<endl;	
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
				if(kappa[ix][iL]>=1.0)
					omega_b_3D[ix][iL][iE] = M_PI *sqrt(kappa[ix][iL])/K[ix][iL]*sqrt(x*tok->eps*L*0.5)/q *sqrt(E);//kappa>1
				else
					omega_b_3D[ix][iL][iE] = M_PI /K[ix][iL]*sqrt(x*tok->eps*L*0.5)/q *sqrt(E);//kappa<1
					
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
//		omega_phi_3D[ix][iL][iE] = -1.0*q_1D[ix] * omega_b_3D[ix][iL][iE];
		omega_phi_3D[ix][iL][iE] = 1.0*(q_1D[ix]) * omega_b_3D[ix][iL][iE];
//		omega_phi_3D[ix][iL][iE] = sqrt(grid->Earray[iE]);
	
	}	

/* for(int ix=0;ix<grid->nx;ix++)
{
        cout<<"omega_b \t"<<omega_b_3D[ix][10][grid->nE-1]<<"\tomega_phi\t"<<omega_phi_3D[ix][10][grid->nE-1]<<endl;
}*/

}




