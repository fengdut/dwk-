#include"math.h"
#include"tokamak.h"
#include"omp.h"
#include<iostream>


void F0_3D(const Slowing* slow,const Grid *grid, double *** F_3D)
{
int np=4;

using namespace std;
	omp_set_num_threads(np);


#pragma omp parallel firstprivate (slow, grid)
{
	int id=omp_get_thread_num();
	
	cout<<"id:\t"<<id<<endl;
//	#pragma omp for
	for(int ix=0;ix<grid->nx;ix=ix++)	
	{
		int nL,nE;
		nL=grid->nL;
		nE=grid->nE;	
		double Ec32=pow(slow->Ec,1.5);
		double x;
		x=grid->dr *ix +grid->ra;
	for(int iL=id;iL<nL;iL=iL+np)
	{
		double L;
		L=grid->dL *iL +grid->La;
	for(int iE=0;iE<nE;iE++)
	{
		double E;
		double TF=0;
		E=grid->dE *iE +grid->Ea;
		
		TF = erfc((E-slow->E0)/slow->Ed)/(pow(E,1.5)+Ec32);
//		double tL=(L-slow->L0)/slow->Ld;
//		tL=-1.0*tL*tL;
//		double tx=(x -slow->r0)/slow->rd;
//		tx=-1.0*tx*tx;
//		F_3D[ix][iL][iE] =  TF*exp(tL)*exp(tx); 
				
	}	
	}
	}
}
}

void F0_3D_test(const Slowing* slow,const Grid *grid, double *** F_3D)
{

	int np=4;
        omp_set_num_threads(np);
	
	int n=grid->nx*grid->nL*grid->nE;

#pragma omp parallel firstprivate(slow,grid) 
{  
	int id=omp_get_thread_num();	
        for(int i=id;i<n;i=i+np)  
        {
                double Ec32=pow(slow->Ec,1.5);
                double x;
                double L;
                double E;
                double TF=0;
		int ix,iL,iE;
		ix=i/(grid->nL*grid->nE);
		iL=(i-ix*grid->nL*grid->nE)/grid->nE;
		iE=i -ix*grid->nL*grid->nE - iL*grid->nE;

               // x=grid->xarray[ix];
		x=grid->dr *ix +grid->ra;
              //  L=grid->Larray[iL];
		L=grid->dL *iL +grid->La;
              //  E=grid->Earray[iE];
		E=grid->dE *iE +grid->Ea;
                
                TF = erfc((E-slow->E0)/slow->Ed)/(pow(E,1.5)+Ec32);
                double tL=(L-slow->L0)/slow->Ld;
                tL=-1.0*tL*tL;
                double tx=(x -slow->r0)/slow->rd;
                tx=-1.0*tx*tx;
                //F_3D[ix][iL][iE] =  TF*exp(tL)*exp(tx);               
		double F3D=0;
                F3D =  TF*exp(tL)*exp(tx);               
                                
        }
}
}

 

