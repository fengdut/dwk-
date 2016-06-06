#include"math.h"
#include"tokamak.h"
#include"omp.h"
#include<iostream>
#include"AllocArray.h"
#include"simpintegral.h"
#include"vector.h"
#include"outlog.h"

using namespace std;
void F0_3D(const Slowing* slow,const Grid *grid,Tokamak *tok, double const rho_h,double const *q_1D,int m,
		double *** F_3D, double ***FE_3D, double ***FR_3D,double *** omega_star,double * Cbeta)
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
		if(slow->rflag==0)
		{
			//exponent profile
			double tx=(x -slow->r0)/slow->rd;
                	tx=-1.0*tx*tx;	
			expx[ix] =exp(tx);
			expx1[ix]=2*expx[ix]*(slow->r0 -x)/(slow->rd*slow->rd);
		}
		else if(slow->rflag==1)
		{
			//polynomin profile
			const double * pc=slow->rc;
			double x2=x*x;
			double x3=x*x2;
			double x4=x2*x2;
			double x5=x4*x;
			double x6=x3*x3;
			double x7=x6*x;
			double x8=x4*x4;
			expx[ix] = pc[0] + pc[1]*x +pc[2]*x2 +pc[3]*x3 +pc[4]*x4 +pc[5]*x5 +pc[6]*x6 +pc[7]*x7 +pc[8]*x8;
			expx1[ix]= pc[1] +2*pc[2]*x+3*pc[3]*x2+4*pc[4]*x3+5*pc[5]*x4+6*pc[6]*x5+7*pc[7]*x6+8*pc[8]*x7; 	
		}
		else
		{
			cerr<<"this option is not available"<<endl;
			exit(1);
		}

		
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
	double **EL;
	Alloc2D(EL,grid->nL,grid->nE);
	for(int iL=0;iL<grid->nL;iL++)
        for(int iE=0;iE<grid->nE;iE++)
	{
		double L,E;
		E=grid->dE *iE +grid->Ea;
                L=grid->dL *iL +grid->La;
		EL[iL][iE] = 2*L*(L - slow->L0)/((slow->Ld*slow->Ld)*E); 
	//	EL[iL][iE] =1;
		
	}

	
	for(int ix=0;ix<grid->nx;ix++) 
	for(int iL=0;iL<grid->nL;iL++)
	for(int iE=0;iE<grid->nE;iE++)
	{
		F_3D[ix][iL][iE]  =expx[ix]*expL[iL]*erfcE[iE];
		FE_3D[ix][iL][iE] =-1.0*expx[ix]*expL[iL]*(erfcE1[iE]+expE[iE]-EL[iL][iE]*erfcE[iE]);
		//FE_3D[ix][iL][iE] =-1.0*expx[ix]*expL[iL]*(erfcE1[iE]);
		//FE_3D[ix][iL][iE] =-1.0*expx[ix]*expL[iL]*(expE[iE]);
		//FE_3D[ix][iL][iE] =-1.0*expx[ix]*expL[iL]*(-1.0*EL[iL][iE]*erfcE[iE]);
		FR_3D[ix][iL][iE]  =expx1[ix]*expL[iL]*erfcE[iE];
		//omega_star[ix][iL][iE] = FR_3D[ix][iL][iE]*tok->R0*rho_h/tok->a *m  /(2.0*grid->xarray[ix]); 
		omega_star[ix][iL][iE] = FR_3D[ix][iL][iE]*tok->R0*rho_h/tok->a *m  /(2.0*grid->xarray[ix]) *  q_1D[ix]; 
		
		if((FE_3D[ix][iL][iE])!=0.0)	
			omega_star[ix][iL][iE] = omega_star[ix][iL][iE]/FE_3D[ix][iL][iE];
		else
			omega_star[ix][iL][iE] =0;
		
	}

	double Cn=0,CP=0;
//	cout<<"omega_star: ";
//	max_min_3D(grid->nx, grid->nL, grid->nE,omega_star);
	//cout<<"F_3D: ";
 	//max_min_3D(grid->nx, grid->nL, grid->nE,F_3D);
	//cout<<"FE_3D: ";
	//max_min_3D(grid->nx, grid->nL, grid->nE,FE_3D);
	//cout<<"expx: ";
	//max_min_1D(grid->nx,expx);
	//cout<<"expL: ";
	//max_min_1D(grid->nL,expL);
	//cout<<"expE: ";
	//max_min_1D(grid->nE,erfcE);
	double **P0,**Th;
	Alloc2D(P0,grid->nL,grid->nE);
	Alloc2D(Th,grid->nL,grid->nE);
        for(int iL=0;iL<grid->nL;iL++)
        for(int iE=0;iE<grid->nE;iE++)
	{
		P0[iL][iE]=F_3D[0][iL][iE]*sqrt(grid->Earray[iE])*sqrt(2)*M_PI/(sqrt(1-grid->Larray[iL]));
		Th[iL][iE]=F_3D[0][iL][iE]*sqrt(grid->Earray[iE])*sqrt(2)*M_PI/(sqrt(1-grid->Larray[iL])) *2*grid->Earray[iE];
	}
	Cn=simpintegral_2D(P0, grid->nL, grid->dL, grid->nE,grid->dE);
	CP=simpintegral_2D(Th, grid->nL, grid->dL, grid->nE,grid->dE);
	
	cout<<"Cn= "<<Cn<<std::endl;
	cout<<"CP= "<<CP<<std::endl;
	*Cbeta = CP/Cn;
	cout<<"Cbeta= "<<*Cbeta<<std::endl;
	cout<<"Th= "<<*Cbeta <<" code unit (Ei_0)"<<std::endl;
	cout<<"E0: "<<tok->E_i0<<std::endl;;
	cout<<"Th= "<<*Cbeta*tok->E_i0 <<" KeV"<<std::endl;
	
	double T2=2*sqrt(2);
	for(int ix=0;ix<grid->nx;ix++)
        for(int iL=0;iL<grid->nL;iL++)
        for(int iE=0;iE<grid->nE;iE++)
	{
		F_3D[ix][iL][iE]  =F_3D[ix][iL][iE]/Cn*T2;
                FE_3D[ix][iL][iE] =FE_3D[ix][iL][iE] /Cn*T2;
                FR_3D[ix][iL][iE]  =FR_3D[ix][iL][iE]/Cn*T2;
	}
	Free2D(P0);
	Free2D(Th);
	Free2D(EL);
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
		//Theta_3D[ix][iL][it] = Theta_3D[ix][iL][it-1] +(1 +1)*grid->dtheta*0.5;
	}
//	cout<<"Theta_3D";
	//for(int ix=0;ix<grid->nx;ix++)
        //for(int iL=0;iL<grid->nL;iL++)
        //for(int it=2;it<grid->ntheta;it++)
//	{
		//Theta_3D[ix][iL][it] = Theta_3D[ix][iL][it-1] +Theta_3D[ix][iL][it];
//	}	
//	max_min_3D(grid->nx,grid->nL,grid->ntheta,Theta_3D);
}





