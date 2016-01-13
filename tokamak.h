#ifndef TOKAMAK_H
#define TOKAMAK_H


int const nqc =5;
typedef struct TOKAMAK
{
	double a,R0;	//minor and major radius (m). 			input
	double eps;	//a/R0	  
	double C; 	//the normalization coefficient. useless now. 	input
	double qc[nqc];   //q profile, polynomial coefficient;      	input
	double s;	//magnetic shear at resonance surface.
	double Bt;	//toroidal magnetic field at axis without plasma, in units of tesla  input
	double Bps;	//poloidal magnetic field at resonance surface.
	double n0;	//thermal plasma density at axis.		input
	double mi;	//ion mass in units of protom mass.		input
	double E_i0;	// fast ion injection energy, in units of KeV	input
	double v_A,omega_A; 
	double tau_At; // tau_{A,theta}
	double v_i0,omega_i0;
	
	
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
	double rho_d;
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

void qprofile(int const nx,double const *xarray, double *q_1D);
void qprofile(int const nx,double const *xarray,double const  qc[5], double *q_1D);

double bf(const Tokamak *tok,const double theta, const double r);
void J_q(Grid *const grid,double *const q_1D, double *J_q);



#endif


