#include"tokamak.h"
#include"math.h"

#include<iostream>
#include<fstream>
#include<ctime>
#include <libconfig.h++>
#include"readinput.h"
#include"simpintegral.h"
#include"AllocArray.h"
#include"vector.h"
#include"dist.h"
#include"mode.h"
#include"particle_fre.h"
#include<complex>
#include"dwk.h"

using namespace std;
void help()
{
	cout<<"test delta wk"<<endl;
}

int main()
{	
	help();
	Tokamak tok;
	Grid 	grid;
	Slowing slowing;

	char filename[]="dwk.cfg";
//read input parameters
	read_tokamak(filename,&tok,&grid,&slowing);	

 	cout.precision (12);
	cout<<"a: \t"<<tok.a<<"\tR0: \t"<<tok.R0<<"\teps: \t"<<tok.eps<<endl;
	cout<<"r0, rd: \t"<<slowing.r0<<"\t"<<slowing.rd<<endl;
	cout<<"L0, Ld: \t"<<slowing.L0<<"\t"<<slowing.Ld<<endl;
	cout<<"E0, Ed, Ec: \t"<<slowing.E0<<"\t"<<slowing.Ed<<"\t"<<slowing.Ec<<endl;	

	cout<<"nx, nL, nE, ntheta: \t"<<grid.nx<<"\t"<<grid.nL<<"\t"<<grid.nE<<"\t"<<grid.ntheta<<endl;
	cout<<"ra, rb: \t"<<grid.ra<<"\t"<<grid.rb<<endl;
	cout<<"La, Lb: \t"<<grid.La<<"\t"<<grid.Lb<<endl;
	cout<<"Ea, Eb: \t"<<grid.Ea<<"\t"<<grid.Eb<<endl;
	
//	cout<<"pi: \t"<<M_PI<<endl;
	clock_t c_start=clock();
	double **test;
Alloc2D(test,grid.nx,grid.nL);
Free2D(test);
	
//to begin to calculate non-omega parts
	double *xarray,*Larray,*Earray, *tarray;
	Alloc1D(xarray,grid.nx);
	Alloc1D(Larray,grid.nL);
	Alloc1D(Earray,grid.nE);
	Alloc1D(tarray,grid.ntheta);
	linspace(grid.ra,grid.rb,grid.nx,xarray);	
	linspace(grid.La,grid.Lb,grid.nL,Larray);
	linspace(grid.Ea,grid.Eb,grid.nE,Earray);
	linspace(0,2*M_PI,grid.ntheta,tarray);
	grid.xarray=xarray;
	grid.Larray=Larray;
	grid.Earray=Earray;
	grid.tarray=tarray;
	
	double *q_1D;
	Alloc1D(q_1D,grid.nx);
	qprofile(grid.nx,xarray,q_1D);
	//ofstream fout("test.dat");
	//printv(fout,grid.nx,q_1D);	
	max_min_1D(grid.nx,q_1D);
	
	double ***F_3D,***FE_3D,***FR_3D, ***omega_star_3D;
	Alloc3D(F_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FE_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FR_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_star_3D,grid.nx,grid.nL,grid.nE);
	double rho_h =0.1;// orbit width
	double m=1;
		
	
	F0_3D(&slowing,&grid,&tok,rho_h,m,F_3D,FE_3D,FR_3D,omega_star_3D);	
	cout<<"F0"<<endl;
	max_min_3D(grid.nx,grid.nL,grid.nE,F_3D);
	cout<<"FE"<<endl;
	max_min_3D(grid.nx,grid.nL,grid.nE,FE_3D);
	cout<<"FR"<<endl;
	max_min_3D(grid.nx,grid.nL,grid.nE,FR_3D);
	cout<<"omega_star_3D: ";
	max_min_3D(grid.nx,grid.nL,grid.nE,omega_star_3D);

	double *** lambda_b_3D, ***b_lambda_3D, ***Theta_3D;
	Alloc3D(lambda_b_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(b_lambda_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(Theta_3D,grid.nx,grid.nL,grid.ntheta);

	Lambda_b_L_3D(&grid,&tok,lambda_b_3D,b_lambda_3D);

        cout<<"lambda_b_3D"<<endl;
        max_min_3D(grid.nx,grid.nL,grid.ntheta,lambda_b_3D);
	
	Theta(b_lambda_3D,&grid,Theta_3D);
	cout<<"dtheta:\t"<<grid.dtheta<<endl;
	cout<<"Theta_3D"<<endl;
	max_min_3D(grid.nx,grid.nL,grid.ntheta,Theta_3D);	

	cout<<"b_lambda_3D"<<endl;
        max_min_3D(grid.nx,grid.nL,grid.ntheta,b_lambda_3D);
	
	double **G_2D;
	Alloc2D(G_2D,grid.nx,grid.ntheta);
	G_R_theta(&grid,&tok,0.5,0.06,G_2D); //r_s=0.5, delta_r=0.06
	cout<<"G_2D:\t";
	max_min_2D(grid.nx,grid.ntheta,G_2D);	
	
	double **Chi_2D, **K_2D,**kappa_2D;
	Alloc2D(Chi_2D,grid.nx,grid.nL);
	Alloc2D(K_2D,grid.nx,grid.nL);
	Alloc2D(kappa_2D,grid.nx,grid.nL);
	
	double sigma=1.0;
	Chi(&grid,&tok,sigma,Chi_2D,kappa_2D,K_2D);	
	cout<<"Chi_2D:";
	max_min_2D(grid.nx,grid.nL,Chi_2D);
		
	cout<<"kappa_2D";
	max_min_2D(grid.nx,grid.nL,kappa_2D);
	
	cout<<"K_2D";
	max_min_2D(grid.nx,grid.nL,K_2D);
	
	complex<double> **Yps_2D;	
	double **Yp2_2D;	
	Alloc2D(Yps_2D,grid.nx,grid.nL);
	Alloc2D(Yp2_2D,grid.nx,grid.nL);
	int p=0; ////////////////////////////////////p number
	Yps(&grid,G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,p,Yps_2D);
	cout<<"Yps_2D:";
	max_min_2D(grid.nx,grid.nL,Yps_2D);
	Yp_2(Yps_2D,&grid,Yp2_2D);
	cout<<"Yp2_2D: ";
	max_min_2D(grid.nx,grid.nL,Yp2_2D);	

	double *** omega_b_3D,***omega_phi_3D,*** tau_b_3D;
	Alloc3D(omega_b_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_phi_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(tau_b_3D,grid.nx,grid.nL,grid.nE);
	omega_b(&grid, &tok,kappa_2D, K_2D, q_1D,omega_b_3D);
	cout<<"omega_b: ";
	max_min_3D(grid.nx,grid.nL,grid.nE,omega_b_3D);

	omega_phi(&grid,q_1D,omega_b_3D,omega_phi_3D);
	cout<<"omega_phi: ";
	max_min_3D(grid.nx,grid.nL,grid.nE,omega_phi_3D);

	tau_b(&grid,omega_b_3D,tau_b_3D);
	cout<<"tau_b: ";
	max_min_3D(grid.nx,grid.nL,grid.nE,tau_b_3D);
	double *J_q_1D;
	Alloc1D(J_q_1D,grid.nx);
	J_q(&grid,q_1D,J_q_1D);
	cout<<"J_q_1D: ";
	max_min_1D(grid.nx,J_q_1D);
	
	complex<double> ***dwk_3D;
	Alloc3D(dwk_3D,grid.nx,grid.nL,grid.nE);
	
	complex<double> dwk_sum;
	complex<double> test_omega=1.5 +0.001i;
	dwk_sum = dwk(&grid, test_omega,  1, 0, omega_phi_3D, omega_b_3D, tau_b_3D, Yp2_2D, J_q_1D, FE_3D, omega_star_3D,dwk_3D);
	cout<<"dwk_3D";
	max_min_3D(grid.nx,grid.nL,grid.nE,dwk_3D);	
	cout<<"dwk(omega): "<<dwk_sum<<endl;	
	
	
	Free3D(dwk_3D);
	Free1D(J_q_1D);
	Free3D(tau_b_3D);
	Free3D(omega_phi_3D);
	Free3D(omega_b_3D);
	Free2D(Yp2_2D);
	Free2D(Yps_2D);
	Free2D(kappa_2D);
	Free2D(K_2D);
	Free2D(Chi_2D);	
	Free2D(G_2D);
	Free3D(Theta_3D);
	Free3D(b_lambda_3D);
	Free3D(lambda_b_3D);
	Free3D(F_3D);	
	Free3D(FR_3D);	
	Free3D(FE_3D);	
	Free3D(omega_star_3D);
	Free1D(q_1D);

	Free1D(tarray);
	Free1D(xarray);
	Free1D(Larray);
	Free1D(Earray);
	
	cout<<"time:"<<float(clock() - c_start)/CLOCKS_PER_SEC<<endl;
}


