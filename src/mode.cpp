#include<iostream>
#include<cmath>
#include"tokamak.h"
#include"AllocArray.h"
#include"vector.h"
#include"mode.h"
using namespace std;

double fxi_r(double x,double r_s, double delta_r)
{
	if(x<=r_s-delta_r*0.5)
	{	
	
		return 1.0;
	}
	if(x>=r_s+delta_r*0.5)
	{
		return 0.0;
	}
	double f=(delta_r -  x  +r_s - 0.5*delta_r)/delta_r;
	return  f;
}
double fxi_t(double x,double r_s, double delta_r)
{
	if(x<=r_s-delta_r*0.5)
	{
                return -1.0*(x);
	}
        if(x>=r_s+delta_r*0.5)
	{
                return 0.0;
	}
	return -x*(delta_r -2*x +r_s - 0.5*delta_r)/delta_r;
}
void G_R_theta(Grid * const grid, Tokamak * const tok, Slowing *const slow,Mode * const pmode, double * const q_1D, 
	complex<double> ***G_3D)
{	
	double r_s=tok->r_s;
	double delta_r =pmode->delta_r;
	double *cost,*sint;
	complex<double> *expt;
	double *xi_r,*xi_t;
	Alloc1D(cost,grid->ntheta);
	Alloc1D(sint,grid->ntheta);
	Alloc1D(expt,grid->ntheta);
	Alloc1D(xi_r,grid->nx);
	Alloc1D(xi_t,grid->nx);
	cout<<"end alloc1d"<<endl;
	
	complex<double> ti=-1.0i;
	for(int it=0;it<grid->ntheta;it++)
	{
		cost[it]=cos(it*grid->dtheta);
		sint[it]=sin(it*grid->dtheta);
		expt[it]=exp(ti*double(it)*grid->dtheta);
	}
	int nx1,nx2;
	nx1 = ceil((r_s-0.5*delta_r)/grid->dr);
	nx2 = ceil((r_s+0.5*delta_r)/grid->dr);
	for(int ix=0;ix<nx1;ix++)
	{
		xi_r[ix]=1.0;
		xi_t[ix]=-1.0;
	}
	for(int ix=nx1;ix<nx2;ix++)
	{
		xi_r[ix] =  (delta_r -  grid->xarray[ix] +r_s - 0.5*delta_r)/delta_r;
		xi_t[ix] = -(delta_r -2*grid->xarray[ix] +r_s - 0.5*delta_r)/delta_r;	
	}		
	for(int ix=nx2;ix<grid->nx;ix++)
	{
		xi_r[ix]=0;
		xi_t[ix]=0;
	}
	for(int ix=0;ix<grid->nx;ix++)
		xi_t[ix] = xi_t[ix] *grid->xarray[ix];
	
	for(int ix=0;ix<grid->nx;ix++)
	{	
		for(int iE=0;iE<grid->nE;iE++)
		{
		//	double rho_d= q_1D[ix] *slow->rho_h*sqrt(grid->Earray[iE]);   //Lambda_0==0 
		for(int it=0;it<grid->ntheta;it++)
		{
			double ttheta   =grid->tarray[it];
			double tgrr,tgtt,tgrt,tkappa_r,tkappa_t;
			complex<double> txi_t,txi_r;
			double b = 1 +	grid->xarray[ix]*tok->eps *cost[it];
			double lb = slow->L0/b;
			double rho_d = q_1D[ix] *0.5 *slow->rho_h *sqrt(grid->Earray[iE]/(1-lb))*(2-lb);
			double x	=grid->xarray[ix]+rho_d*cost[it];	
			double teps	=tok->eps*x;
			
			tgrr =1.0 +teps *cost[it]*0.5;
			tgtt =(1.0-2.5 *teps *cost[it])/(x*x);
			tgrt =-1.5*teps *sint[it] /x;
			
			tkappa_r = -1.0*cost[it]*tok->eps +tok->eps *teps*0.25 - tok->eps *1.25*teps*(cos(2*ttheta)-1)-tok->eps*tok->eps *x/q_1D[ix];
			tkappa_t = teps*sint[it]+ 1.25 *teps*teps *sin(2*ttheta);

			txi_r  = 	 expt[it]*fxi_r(grid->xarray[ix]+rho_d*cost[it],r_s,delta_r);
			txi_t  = -1.0*ti*expt[it]*fxi_t(grid->xarray[ix]+rho_d*cost[it],r_s,delta_r);
			
			G_3D[ix][iE][it] =    (tgtt *tkappa_t +tgrt *tkappa_r)*txi_t + (tgrr *tkappa_r +tgrt *tkappa_t)*txi_r;
					
		}	
		}
	}
/*	
	cout<<xi_r2d[0][4]<<endl;
	cout<<"grr:\t";
	max_min_2D(grid->nx,grid->ntheta,grr);
	cout<<"gtt:\t";
	max_min_2D(grid->nx,grid->ntheta,gtt);
	cout<<"grt:\t";
	max_min_2D(grid->nx,grid->ntheta,grt);

     	cout<<"kappa_r:\t";
        max_min_2D(grid->nx,grid->ntheta,kappa_r);
	
     	cout<<"kappa_t:\t";
        max_min_2D(grid->nx,grid->ntheta,kappa_t);

	cout<<"xi_t:\t";
        max_min_2D(grid->nx,grid->ntheta,xi_t2d);

	cout<<"xi_r:\t";
        max_min_2D(grid->nx,grid->ntheta,xi_r2d);
*/	
	Free1D(xi_t);
	Free1D(xi_r);
	Free1D(cost);
	Free1D(sint);	
	Free1D(expt);
}





