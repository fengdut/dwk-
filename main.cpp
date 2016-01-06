#include<fstream>
#include<iostream>
#include"math.h"
#include<ctime>
#include<complex>
#include<libconfig.h++>

#include"vector.h"
#include"AllocArray.h"
#include"tokamak.h"
#include"readinput.h"
#include"dist.h"
#include"mode.h"
#include"particle_fre.h"
#include"dwk.h"

using namespace std;
int main()
{	
	clock_t c_start=clock();
	help();
	Tokamak tok;
	Grid 	grid;
	Slowing slowing;
	Mode mode;
	char filename[]="dwk.cfg";
//read input parameters
	read_tokamak(filename,&tok,&grid,&slowing,&mode);	
	CGrid pgrid(&grid,&slowing);
 	cout.precision (12);
	showtokamak(&tok,&slowing);
	pgrid.showgrid();	
	
	double *q_1D,*J_q_1D;
	Alloc1D(q_1D,grid.nx);
	Alloc1D(J_q_1D,grid.nx);
	double **G_2D, **Chi_2D, **K_2D,**kappa_2D, **Yp2_2D;	
	Alloc2D(G_2D,grid.nx,grid.ntheta);
	Alloc2D(Chi_2D,grid.nx,grid.nL);
	Alloc2D(K_2D,grid.nx,grid.nL);
	Alloc2D(kappa_2D,grid.nx,grid.nL);
	Alloc2D(Yp2_2D,grid.nx,grid.nL);
	complex<double> **Yps_2D;	
	Alloc2D(Yps_2D,grid.nx,grid.nL);
	double *** omega_b_3D,***omega_phi_3D,*** tau_b_3D;
	Alloc3D(omega_b_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_phi_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(tau_b_3D,grid.nx,grid.nL,grid.nE);
	double ***F_3D,***FE_3D,***FR_3D, ***omega_star_3D;
	Alloc3D(F_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FE_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FR_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_star_3D,grid.nx,grid.nL,grid.nE);
	double *** lambda_b_3D, ***b_lambda_3D, ***Theta_3D;
	Alloc3D(lambda_b_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(b_lambda_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(Theta_3D,grid.nx,grid.nL,grid.ntheta);
	complex<double> ***dwk_3D;
	Alloc3D(dwk_3D,grid.nx,grid.nL,grid.nE);

//begin calculate non-omega parts------------------------
	qprofile(grid.nx,grid.xarray,q_1D);
	F0_3D(&slowing,&grid,&tok,slowing.rho_h,mode.m,F_3D,FE_3D,FR_3D,omega_star_3D);	
	Lambda_b_L_3D(&grid,&tok,lambda_b_3D,b_lambda_3D);
	Theta(b_lambda_3D,&grid,Theta_3D);
	G_R_theta(&grid,&tok,&slowing,&mode,G_2D); 
	Chi(&grid,&tok,slowing.sigma,Chi_2D,kappa_2D,K_2D);	
	omega_b(&grid, &tok,kappa_2D, K_2D, q_1D,omega_b_3D);
	omega_phi(&grid,q_1D,omega_b_3D,omega_phi_3D);
	tau_b(&grid,omega_b_3D,tau_b_3D);
	J_q(&grid,q_1D,J_q_1D);
//end calculate non-omega parts--------------------------	
	complex<double> dwk_sum=0;
	ofstream fout("dwk.out");
 	fout.precision (12);

 	complex<double> tdwk_sum;
        complex<double> err=1;
        int maxi=30;
        complex<double> omega_0=0.0+0.001i;
        complex<double>ti=1i;
        for(int tt=0;tt<maxi;tt++)
	{
		dwk_sum=0;
		tdwk_sum=0;
		omega_0 =omega_0+0.05;
		//omega_0=0;

		for(int ip=mode.pa;ip<=mode.pb;ip++)
		{
			Yps(&grid,G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,ip,Yps_2D);
			Yp_2(Yps_2D,&grid,Yp2_2D);
			tdwk_sum = dwk(&grid, omega_0,mode.n,ip, omega_phi_3D, omega_b_3D, tau_b_3D, Yp2_2D, J_q_1D, FE_3D, omega_star_3D,dwk_3D);
			cout<<"p="<<ip<<"\t dwk_sum:"<<tdwk_sum<<endl;
			dwk_sum =dwk_sum +tdwk_sum;
		}	
	        err=ti*dwk_sum*tok.C +omega_0;
                cout<<"omega_0:"<<omega_0<<"\t dwk_sum:"<<dwk_sum<<"\t err:"<<err<<endl;
	//	omega_0=0.05+imag(omega_0)*ti;

		fout<<real(omega_0)<<"\t"<<imag(omega_0)<<"\t"<<real(dwk_sum)<<"\t"<<imag(dwk_sum)<<endl;	
   //             omega_0=real(-dwk_sum*ti*tok.C)+0.002i;
   //             omega_0=-dwk_sum*ti*tok.C;
//		cout<<"omega_0:"<<omega_0<<endl;
                if(abs(err)<1e-10)
                         break;

	}
	
	Free1D(J_q_1D);		Free1D(q_1D);
        Free2D(Chi_2D); 	Free2D(G_2D);
        Free2D(kappa_2D);	Free2D(K_2D);
	Free2D(Yp2_2D);		Free2D(Yps_2D);
	Free3D(dwk_3D); 	Free3D(tau_b_3D);
	Free3D(omega_phi_3D);	Free3D(omega_b_3D);
	Free3D(Theta_3D);
	Free3D(b_lambda_3D);	Free3D(lambda_b_3D);
	Free3D(F_3D);		Free3D(FR_3D);	
	Free3D(FE_3D);		Free3D(omega_star_3D);
	cout<<"time:\t"<<float(clock() - c_start)/CLOCKS_PER_SEC<<endl;
}


