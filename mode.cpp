
#include<iostream>
#include<cmath>
#include"tokamak.h"
#include"AllocArray.h"
#include"vector.h"
#include"mode.h"
using namespace std;

void G_R_theta(Grid * const grid, Tokamak * const tok,Mode * const pmode,double **G_2D)
{	
	double r_s=pmode->r_s;
	double delta_r =pmode->delta_r;
	double *cost,*sint;
	double *xi_r,*xi_t;
	Alloc1D(cost,grid->ntheta);
	Alloc1D(sint,grid->ntheta);
	Alloc1D(xi_r,grid->nx);
	Alloc1D(xi_t,grid->nx);
	
	for(int it=0;it<grid->ntheta;it++)
	{
		cost[it]=cos(it*grid->dtheta);
		sint[it]=sin(it*grid->dtheta);
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
	
	double ** xi_r2d, ** xi_t2d;
	Alloc2D(xi_r2d,grid->nx,grid->ntheta);
	Alloc2D(xi_t2d,grid->nx,grid->ntheta);
		
	for(int ix=0;ix<grid->nx;ix++)
	for(int it=0;it<grid->ntheta;it++)
	{
		xi_r2d[ix][it] = cost[it]*xi_r[ix];
		xi_t2d[ix][it] = sint[it]*xi_t[ix];
	}
	
	double **grr,**gtt,**grt,**kappa_r,**kappa_t;
	Alloc2D(grr,grid->nx,grid->ntheta);
	Alloc2D(gtt,grid->nx,grid->ntheta);
	Alloc2D(grt,grid->nx,grid->ntheta);
	Alloc2D(kappa_r,grid->nx,grid->ntheta);
	Alloc2D(kappa_t,grid->nx,grid->ntheta);
			
	for(int ix=0;ix<grid->nx;ix++)
	{	
		double x	=grid->xarray[ix];	
		double teps	=tok->eps*grid->xarray[ix];
		for(int it=0;it<grid->ntheta;it++)
		{
			grr[ix][it] =1.0 +teps *cost[it]*0.5;
			gtt[ix][it] =(1.0-2.5 *teps *cost[it])/(x*x);
			grt[ix][it] =-1.5*teps *sint[it] /x;
			kappa_r[ix][it] = -1.0*cost[it]/(1/tok->eps +x *cost[it]);
			kappa_t[ix][it] = teps*sint[it];
			G_2D[ix][it] =    (gtt[ix][it] *kappa_t[ix][it] +grt[ix][it] *kappa_r[ix][it])*xi_t2d[ix][it] 
					+ (grr[ix][it] *kappa_r[ix][it] +grt[ix][it] *kappa_t[ix][it])*xi_r2d[ix][it];
					
		}
	}
	
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
	

	Free2D(kappa_t);
	Free2D(kappa_r);
	Free2D(grt);
	Free2D(gtt);
	Free2D(grr);
	Free2D(xi_r2d);
	Free2D(xi_t2d);
	Free1D(xi_t);
	Free1D(xi_r);
	Free1D(cost);
	Free1D(sint);	
}



