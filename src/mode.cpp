#include<iostream>
#include<cmath>
#include"tokamak.h"
#include"AllocArray.h"
#include"vector.h"
#include"mode.h"
#include"dwk.h"
using namespace std;

double fxi_r(double x,double r_s, double delta_r)
{
	if(x<=r_s-delta_r*0.5)
	//if(x<=r_s)
	{	
	
		return 1.0;
	}
	if(x>=r_s+delta_r*0.5)
	//if(x>=r_s+delta_r)
	{
		return 0.0;
	}
	double f=(delta_r -  x  +r_s - 0.5*delta_r)/delta_r;
	//double f=1.0 -  (x -r_s )/delta_r;
	return  f;
}

double fxi_r_layer0(double x,double r_s, double delta_r)
{
        //if(x<=r_s-delta_r*0.5)
        if(x<=r_s)
        {

                return 1.0;
        }
                return 0.0;
}

double fxi_t(double x,double r_s, double delta_r)
{
	if(x<=r_s-delta_r*0.5)
	//if(x<=r_s)
	{
                return -1.0*(x);
	}
        if(x>=r_s+delta_r*0.5)
        //if(x>=r_s+delta_r)
	{
                return 0.0;
	}
	return -x*(delta_r -2*x +r_s - 0.5*delta_r)/delta_r;
	//return -x*(1 -(2*x -r_s)/delta_r);
//	return -x;
}

double fxi_t_layer0(double x,double r_s, double delta_r)
{
        if(x<=r_s)
        {
                return -1.0*(x);
        }
                return 0.0;
}

void G_R_theta(Grid * const grid, Tokamak * const tok, Slowing *const slow,Mode * const pmode, double * const q_1D, 
	complex<double> ***G_3D)
{	
	double r_s=tok->r_s;
	double delta_r =pmode->delta_r;
	double *cost,*sint;
	complex<double> *expt;
	Alloc1D(cost,grid->ntheta);
	Alloc1D(sint,grid->ntheta);
	Alloc1D(expt,grid->ntheta);
	double (*fxi_t_t)(double, double, double )=fxi_t;
	double (*fxi_r_t)(double, double, double )=fxi_r;
	if(pmode->zero_iner==0)
	{
		(fxi_t_t)=fxi_t_layer0;
		(fxi_r_t)=fxi_r_layer0;
	}
	
	complex<double> ti;
	ti.real(0.0);
	ti.imag(-1.0);
	for(int it=0;it<grid->ntheta;it++)
	{
		cost[it]=cos(it*grid->dtheta);
		sint[it]=sin(it*grid->dtheta);
		expt[it]=exp(ti*double(it)*grid->dtheta);
	}
	int nx1,nx2;
	nx1 = ceil((r_s-0.5*delta_r)/grid->dr);
	nx2 = ceil((r_s+0.5*delta_r)/grid->dr);
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
			//double b = 1 +	grid->xarray[ix]*tok->eps *cost[it];
			double b = 1 +	grid->xarray[ix]*tok->eps;
			double lb = slow->L0/b;
			double rho_d = q_1D[ix] *0.5 *slow->rho_h *sqrt(grid->Earray[iE]/(1-lb))*(2-lb);
	//		double rho_d = slow->rho_h;
			rho_d=rho_d *pmode->zero_rhod;
			double x	=grid->xarray[ix]+rho_d*cost[it];	
			double teps	=tok->eps*x;
			
			tgrr =1.0 +teps *cost[it]*0.5;
			tgtt =(1.0-2.5 *teps *cost[it])/(x*x);
			tgrt =-1.5*teps *sint[it] /x;
			
			tkappa_r = -1.0*cost[it]*tok->eps +tok->eps *teps*0.25 - tok->eps *1.25*teps*(cos(2*ttheta)-1)-tok->eps*tok->eps *x/q_1D[ix];
			tkappa_t = teps*sint[it]+ 1.25 *teps*teps *sin(2*ttheta);

			txi_r  = 	 expt[it]*(*fxi_r_t)(grid->xarray[ix]+rho_d*cost[it],r_s,delta_r)*pmode->xi_0;
			txi_t  = -1.0*ti*expt[it]*(*fxi_t_t)(grid->xarray[ix]+rho_d*cost[it],r_s,delta_r)*pmode->xi_0;
			
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
	Free1D(cost);
	Free1D(sint);	
	Free1D(expt);
}





