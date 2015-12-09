#include"math.h"
#include"tokamak.h"
#include"omp.h"
#include<iostream>
#include"AllocArray.h"

void F0_3D(const Slowing* slow,const Grid *grid, double *** F_3D, double ***FE_3D, double ***FR_3D)
{
	double *expx,*expx1,*expL,*expE, *erfcE,*erfcE1;
	Alloc1D(expx,grid->nx);
	Alloc1D(expx1,grid->nx);
	Alloc1D(expL,grid->nL);
	Alloc1D(expE,grid->nE);	
	Alloc1D(erfcE,grid->nE);
	Alloc1D(erfcE1,grid->nE);
	for(int ix=0;ix<grid->nx;ix++)	
	{
		double x;
		x=grid->dr *ix +grid->ra;
		double tx=(x -slow->r0)/slow->rd;
                tx=-1.0*tx*tx;	
		expx[ix] =exp(tx);
		expx1[ix]=2*expx[ix]*(slow->r0 -x)/(slow->rd*slow->rd);
		
	}	
	for(int iL=0;iL<grid->nL;iL++)
	{
		double L;
		L=grid->dL *iL +grid->La;
		double tL=(L-slow->L0)/slow->Ld;
                tL=-1.0*tL*tL;
		expL[iL] = exp(tL);
	}
	double Ec32=pow(slow->Ec,1.5);
	double sqrtpi =pow(M_PI,0.5);
	for(int iE=0;iE<grid->nE;iE++)
	{
		double E;
		E=grid->dE *iE +grid->Ea;
		double tE = (E-slow->E0)/slow->Ed;
		double E32 =pow(E,1.5)+Ec32;
		erfcE[iE] = erfc(tE)/E32;
		tE=-1.0*tE*tE;
		expE[iE]  = 2*exp(tE)/(sqrtpi*slow->Ed*E32);		
		erfcE1[iE] =1.5*pow(E,0.5)*erfcE[iE]/E32;
		
	}
	for(int ix=0;ix<grid->nx;ix++) 
	for(int iL=0;iL<grid->nL;iL++)
	for(int iE=0;iE<grid->nE;iE++)
	{
		F_3D[ix][iL][iE]  =expx[ix]*expL[iL]*erfcE[iE];
		FE_3D[ix][iL][iE] =-1.0*expx[ix]*expL[iL]*(erfcE1[iE]+expE[iE]);

		FR_3D[ix][iL][iE]  =expx1[ix]*expL[iL]*erfcE[iE];
		
	}

	Free1D(expx);
	Free1D(expx1);
	Free1D(expL);
	Free1D(expE);
	Free1D(erfcE);
	Free1D(erfcE1);
}


void Lambda_b_L_3D(const Grid *grid, const Tokamak *tok,double*** lambda_b_3D,double *** b_lambda_3D)
{
	for(int ix=0;ix<grid->nx;ix++)
	for(int iL=0;iL<grid->nL;iL++)	
	for(int it=0;it<grid->ntheta;it++)
	{
		double L=grid->dL *iL +grid->La;
		double x=grid->dr *ix +grid->ra;
		double theta=grid->dtheta*it;
		double tb =bf(tok,theta,x);
		lambda_b_3D[ix][iL][it] = L/tb;
		b_lambda_3D[ix][iL][it] = 1/(tb *sqrt(1-lambda_b_3D[ix][iL][it]));
	}
}

void Theta(double ***b_lambda_3D,const Grid *grid, double ***Theta_3D)
{
 	for(int ix=0;ix<grid->nx;ix++)
        for(int iL=0;iL<grid->nL;iL++)
        for(int it=1;it<grid->ntheta;it++)
	{
		Theta_3D[ix][iL][it] = Theta_3D[ix][iL][it-1] +(b_lambda_3D[ix][iL][it] +b_lambda_3D[ix][iL][it-1])*grid->dtheta*0.5;
	}
	for(int ix=0;ix<grid->nx;ix++)
        for(int iL=0;iL<grid->nL;iL++)
        for(int it=2;it<grid->ntheta;it++)
	{
		Theta_3D[ix][iL][it] = Theta_3D[ix][iL][it-1] +Theta_3D[ix][iL][it];
	}	
}





