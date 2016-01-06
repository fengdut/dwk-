#include<iostream>
#include<assert.h>

#include"cmath"
#include"vector.h"
#include"AllocArray.h"
#include"tokamak.h"
using namespace std;

void qprofile(const int nx,const double *xarray, double *q_1D)
{
	for(int i=0;i<nx;i++)
		q_1D[i] = 0.9 +2*xarray[i]*xarray[i];
	
}


double bf(const Tokamak *tok,const double theta, const double r)
{
	return 1-tok->eps*r*cos(theta);
}

void J_q(Grid *const grid,double *const q_1D, double *J_q)
{
	for(int ix=0;ix<grid->nx;ix++)
	{
		J_q[ix] = grid->xarray[ix]/q_1D[ix];
	}	
}


CGrid::CGrid()
{
	m_pgrid=0;
}
CGrid::CGrid(Grid *pgrid,Slowing * pslowing )
{
	m_pgrid=pgrid;
	assert(m_pgrid->nx>0&&m_pgrid->nL>0&&m_pgrid->nE>0);

	Alloc1D(xarray,m_pgrid->nx);
        Alloc1D(Larray,m_pgrid->nL);
        Alloc1D(Earray,m_pgrid->nE);
        Alloc1D(tarray,m_pgrid->ntheta);	

	m_pgrid->xarray=xarray;
        m_pgrid->Larray=Larray;
        m_pgrid->Earray=Earray;
        m_pgrid->tarray=tarray;
	
        m_pgrid->ra=1e-15;
        m_pgrid->rb=1;
        m_pgrid->dr= (pgrid->rb- pgrid->ra)/(pgrid->nx-1);

        m_pgrid->La=1e-9;
        m_pgrid->Lb=pslowing->L0+pslowing->Ld+0.1;
        m_pgrid->dL= (pgrid->Lb - pgrid->La)/(pgrid->nL-1);

        m_pgrid->Ea=0.1;
        m_pgrid->Eb=pslowing->E0*1.2;
        m_pgrid->dE= (pgrid->Eb -pgrid->Ea)/(pgrid->nE-1);

        m_pgrid->dtheta=2*M_PI/(pgrid->ntheta-1);
	
	linspace(m_pgrid->ra,m_pgrid->rb,m_pgrid->nx,xarray);
        linspace(m_pgrid->La,m_pgrid->Lb,m_pgrid->nL,Larray);
        linspace(m_pgrid->Ea,m_pgrid->Eb,m_pgrid->nE,Earray);
        linspace(0,2*M_PI,m_pgrid->ntheta,tarray);
	
}

CGrid::~CGrid()
{
	if(m_pgrid!=0)
	{
	        Free1D(tarray);
        	Free1D(xarray);
        	Free1D(Larray);
        	Free1D(Earray);
	}
}


void CGrid::showgrid()
{
	cout<<"--------begin grid information--------"<<endl;
	cout<<"nr, nL, nE, ntheta:\t"<<m_pgrid->nx<<", "<<m_pgrid->nL<<", "<<m_pgrid->nE<<", "<<m_pgrid->ntheta<<endl;
        cout<<"ra, rb: \t"<<m_pgrid->ra<<"\t"<<m_pgrid->rb<<endl;
        cout<<"La, Lb: \t"<<m_pgrid->La<<"\t"<<m_pgrid->Lb<<endl;
        cout<<"Ea, Eb: \t"<<m_pgrid->Ea<<"\t"<<m_pgrid->Eb<<endl;
	cout<<"--------end grid information  --------"<<endl;
	
}
void showtokamak(Tokamak *ptok,Slowing *pslowing)
{
	cout<<"--------begin tokamak information--------"<<endl;
  	cout<<"a: \t"<<ptok->a<<endl;
	cout<<"R0: \t"<<ptok->R0<<endl;
	cout<<"eps: \t"<<ptok->eps<<endl;
	cout<<"--------end tokamak information  --------"<<endl;	

	cout<<"--------begin fast ion distribution function--------"<<endl;
        cout<<"r0, rd: \t"<<pslowing->r0<<"\t"<<pslowing->rd<<endl;
        cout<<"L0, Ld: \t"<<pslowing->L0<<"\t"<<pslowing->Ld<<endl;
        cout<<"E0, Ed, Ec: \t"<<pslowing->E0<<"\t"<<pslowing->Ed<<"\t"<<pslowing->Ec<<endl;
	cout<<"--------end fast ion distribution function  --------"<<endl;
		
}

