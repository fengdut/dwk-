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
	complex<double> **G_2D)
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
	
	complex<double> ** xi_r2d, ** xi_t2d;
	Alloc2D(xi_r2d,grid->nx,grid->ntheta);
	Alloc2D(xi_t2d,grid->nx,grid->ntheta);
		
	complex<double> ti=-1.0i;
	for(int ix=0;ix<grid->nx;ix++)
	for(int it=0;it<grid->ntheta;it++)
	{
	//	cout<<"fxi_r:"<<grid->xarray[ix]+slow->rho_d*cost[it]<<"\t"<<fxi_r(grid->xarray[ix]+slow->rho_d*cost[it],r_s,delta_r)<<endl;
		xi_r2d[ix][it] = exp(ti*(double)it*grid->dtheta)   *fxi_r(grid->xarray[ix]+slow->rho_d*cost[it],r_s,delta_r);  //rho_d should be here
		xi_t2d[ix][it] = -1.0*ti*exp(ti*(double)it*grid->dtheta)*fxi_t(grid->xarray[ix]+slow->rho_d*cost[it],r_s,delta_r);  //rho_d should be here
	}
	cout<<"?????"<<endl;
	cout<<xi_r2d[0][4]<<endl;
	
//	cout<<"xi_2D ";
//	max_min_2D(grid->nx,grid->ntheta,xi_r2d);
//	max_min_2D(grid->nx,grid->ntheta,xi_t2d);

	double **grr,**gtt,**grt,**kappa_r,**kappa_t;
	Alloc2D(grr,grid->nx,grid->ntheta);
	Alloc2D(gtt,grid->nx,grid->ntheta);
	Alloc2D(grt,grid->nx,grid->ntheta);
	Alloc2D(kappa_r,grid->nx,grid->ntheta);
	Alloc2D(kappa_t,grid->nx,grid->ntheta);
			
	cout<<xi_r2d[0][4]<<endl;
	for(int ix=0;ix<grid->nx;ix++)
	{	
		for(int it=0;it<grid->ntheta;it++)
		{
			double x	=grid->xarray[ix]+slow->rho_d*cost[it];	
			double teps	=tok->eps*x;
			double ttheta   =grid->tarray[it];
			grr[ix][it] =1.0 +teps *cost[it]*0.5;
			gtt[ix][it] =(1.0-2.5 *teps *cost[it])/(x*x);
			grt[ix][it] =-1.5*teps *sint[it] /x;
			

			kappa_r[ix][it] = -1.0*cost[it]*tok->eps +tok->eps *teps*0.25 - tok->eps *1.25*teps*(cos(2*ttheta)-1)-tok->eps*tok->eps *x/q_1D[ix];
			kappa_t[ix][it] = teps*sint[it]+ 1.25 *teps*teps *sin(2*ttheta);
			G_2D[ix][it] =    (gtt[ix][it] *kappa_t[ix][it] +grt[ix][it] *kappa_r[ix][it])*(xi_t2d[ix][it]) 
					+ (grr[ix][it] *kappa_r[ix][it] +grt[ix][it] *kappa_t[ix][it])*(xi_r2d[ix][it]);
					
		}
	}
	
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

void find_rs(Grid * const grid,double *const q_1D,double const q_s, double *r_s)
{
	assert(q_s>0&&q_s<100);
	for(int i=1;i<grid->nx;i++)
	{
		if((q_1D[i-1]-q_s)*(q_1D[i]-q_s)<=0)
		{
			double dr =grid->dr;
			double dq=abs(q_s-q_1D[i-1])/abs((q_1D[i]-q_1D[i-1]));
			*r_s = 	grid->xarray[i-1] +dr *(1-dq);
			return;	
			
		}
	}
	
	cout<<"can not find r_s with q_s= "<<q_s<<endl;
}





