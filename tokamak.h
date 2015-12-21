#ifndef TOKAMAK_H
#define TOKAMAK_H

typedef struct TOKAMAK
{
	double a,R0;
	double eps;
}Tokamak;

typedef struct GRID
{
	int nx,nL,nE,ntheta;
	double ra,rb,dr;
	double La,Lb,dL;
	double Ea,Eb,dE;
	double dtheta;
	double *xarray,*Larray,*Earray,*tarray;
	
}Grid;

typedef struct SLOWING
{
	double  rd,r0;
	double  L0,Ld;
        double 	E0,Ed,Ec;
	double rho_h;
	int sigma;
}
Slowing;


class CGrid
{
public:
	CGrid();
	CGrid(Grid *pgrid,Slowing *slowing);
	~CGrid();
	Grid * m_pgrid;
	void showgrid();
private:
	double *xarray,*Larray,*Earray,*tarray;
	
};


void showtokamak(Tokamak *ptok,Slowing *pslowing);


void qprofile(const int nx,const double *xarray, double *q_1D);
double bf(const Tokamak *tok,const double theta, const double r);
void J_q(Grid *const grid,double *const q_1D, double *J_q);



#endif


