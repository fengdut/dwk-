#include<iostream>
#include<assert.h>

#include"cmath"
#include"vector.h"
#include"AllocArray.h"
#include"tokamak.h"
#include"stdlib.h"
using namespace std;


void qprofile(Grid *const grid,Tokamak *ptok, double *q_1D)
{
	double *qc =ptok->qc;
	for(int i=0;i<grid->nx;i++)
	{
		double x=grid->xarray[i];
		double x2=x*x;
		double x3=x2*x;
		double x4=x3*x;
		double x5=x4*x;
		double x6=x5*x;
		double x7=x6*x;
		q_1D[i] = qc[0] +qc[1]*x +qc[2]*x2 + qc[3]*x3 +qc[4]*x4+qc[5]*x5+qc[6]*x6+qc[7]*x7;
	}
	
	assert(ptok->q_s>0&&ptok->q_s<100);
	ptok->r_s=-1.0;
        for(int i=1;i<grid->nx;i++)
        {
                if((q_1D[i-1]-ptok->q_s)*(q_1D[i]-ptok->q_s)<=0)
                {
                        double dr =grid->dr;
                        double dq=abs(ptok->q_s-q_1D[i-1])/abs((q_1D[i]-q_1D[i-1]));
                        ptok->r_s =  grid->xarray[i-1] +dr *(1-dq);
			break;
                }
        }

	if(ptok->r_s<0)
	{
        	cerr<<"can not find r_s with q_s= "<<ptok->q_s<<endl;
        	exit(-1);
	}	
	double rs=ptok->r_s;
	ptok->s = qc[1] +2*qc[2] *rs +3*qc[3]*rs*rs +4*qc[4]*rs*rs*rs +5*qc[5] *rs*rs*rs*rs +6*qc[6] *rs*rs*rs*rs*rs +7*qc[7] *rs*rs*rs*rs*rs*rs;
	ptok->s =rs *ptok->s;
	
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
        m_pgrid->Lb=pslowing->L0+pslowing->Ld*2.0;
        m_pgrid->dL= (pgrid->Lb - pgrid->La)/(pgrid->nL-1);

        m_pgrid->Ea=1e-2;
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
void calculate_normalization(Tokamak *ptok, Slowing *pslowing)
{
	ptok->eps=ptok->a/ptok->R0;
        ptok->Bps=ptok->r_s*ptok->a*ptok->Bt/(ptok->R0*ptok->q_s);
        ptok->rho_m = ptok->mi * ptok->n0 *1.6726e-27;
        ptok->tau_At =sqrt(3.0) *ptok->r_s*ptok->a / (ptok->Bps/sqrt(ptok->rho_m*4.0*M_PI*1.0e-7));
        ptok->omega_A = 2.0 /(ptok->tau_At *ptok->s);
        ptok->v_i0 = sqrt(2.0*ptok->E_i0 *1.0e3*1.6022e-19/(ptok->m_ep*1.6726e-27));
        ptok->omega_i0 = ptok->v_i0/ptok->R0;


        ptok->C = ptok->omega_A/(ptok->omega_i0*4.0/M_PI *(ptok->r_s*ptok->a/ptok->R0*0.5)*(ptok->r_s*ptok->a/ptok->R0*0.5));
        ptok->beta_h=0.0;

        pslowing->rho_h=ptok->m_ep*1.6726e-27 *ptok->v_i0/(1.6022e-19*ptok->Bt)/ptok->a;
}

void showtokamak(Tokamak *ptok,Slowing *pslowing)
{
	cout<<"--------begin tokamak information--------"<<endl;
  	cout<<"a: \t"<<ptok->a<<" m"<<endl;
	cout<<"R0: \t"<<ptok->R0<<" m"<<endl;
	cout<<"eps: \t"<<ptok->eps<<endl;
	cout<<"Bt0: \t"<<ptok->Bt<<" Tesla"<<endl;
	cout<<"Bps: \t"<<ptok->Bps<<" Tesla"<<endl;
	cout<<"n0: \t"<<ptok->n0<<" 1/m^3"<<endl;
	cout<<"mi: \t"<<ptok->mi<<" 1/m_p"<<endl;
	double *qc=ptok->qc;
	cout<<"q profile: \t"<<qc[0]<<"+"<<qc[1]<<"*r+"<<qc[2]<<"*r^2+"<<qc[3]<<"*r^3+"<<qc[4]<<"*r^4+"<<qc[5]<<"*r^5+"<<qc[6]<<"*r^6+"<<qc[7]<<"*r^7"<<endl;
	cout<<"q_s: \t"<<ptok->q_s<<endl;
	cout<<"r_s: \t"<<ptok->r_s<<endl;
	cout<<"s : \t"<<ptok->s<<endl;
	cout<<"rho_m:\t"<<ptok->rho_m<<" kg/m^3"<<endl;
	cout<<"tau_At:\t"<<ptok->tau_At<<" s"<<endl;
	cout<<scientific;
	cout<<"omega_A:\t"<<ptok->omega_A<< " rad/s"<<endl;
	cout<<fixed;
	cout<<"E_i0:\t"<<ptok->E_i0<<" kEv"<<endl;
	cout<<scientific;
	cout<<"v_i0:\t"<<ptok->v_i0<<" m/s"<<endl;
	cout<<"omega_i0:\t"<<ptok->omega_i0 <<" rad/s"<<endl;
	cout<<"fast ion gyroradius:\t"<<pslowing->rho_h*ptok->a <<" m"<<endl;
	cout<<"C:\t"<<ptok->C<<endl;
	cout<<fixed;
	cout<<"--------end tokamak information  --------"<<endl;	

	cout<<"--------begin fast ion distribution function--------"<<endl;
        cout<<"r0, rd: \t"<<pslowing->r0<<"\t"<<pslowing->rd<<endl;
        cout<<"L0, Ld: \t"<<pslowing->L0<<"\t"<<pslowing->Ld<<endl;
        cout<<"E0, Ed, Ec: \t"<<pslowing->E0<<"\t"<<pslowing->Ed<<"\t"<<pslowing->Ec<<endl;
	cout<<"--------end fast ion distribution function  --------"<<endl;
		
}

void find_rs(Grid * const grid,double *const q_1D, Tokamak *ptok)
{
        assert(ptok->q_s>0&&ptok->q_s<100);
        for(int i=1;i<grid->nx;i++)
        {
                if((q_1D[i-1]-ptok->q_s)*(q_1D[i]-ptok->q_s)<=0)
                {
                        double dr =grid->dr;
                        double dq=abs(ptok->q_s-q_1D[i-1])/abs((q_1D[i]-q_1D[i-1]));
                        ptok->r_s =  grid->xarray[i-1] +dr *(1-dq);
                        return;

                }
        }

        cout<<"can not find r_s with q_s= "<<ptok->q_s<<endl;
        exit(-1);
}

