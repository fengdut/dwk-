#include"tokamak.h"
#include<complex>
#include<cmath>
#include<iostream>
#include"AllocArray.h"
#include"simpintegral.h"
#include"mode.h"

using namespace std;

complex<double> dwk(Grid *const grid, complex<double> const omega, int n,int p, 
	double *** const omega_phi_3D, double ***const omega_b_3D, double *** tau_b_3D,
	double ** const Yp2_2D, double * const J_q_1D, double *** const F_E_3D, double *** const omega_star,complex<double> *** dwk_3D )
{

	double ismall=0.0001;
	if(abs(imag(omega))<ismall)
	{
//		cout<<"image part of omega is too small"<<endl;
	}		
	
	for(int ix=0;ix<grid->nx;ix++)
	for(int iL=0;iL<grid->nL;iL++)
	for(int iE=0;iE<grid->nE;iE++)
	{
		complex<double> Yp_R =Yp2_2D[ix][iL] /(n*omega_phi_3D[ix][iL][iE] +p *omega_b_3D[ix][iL][iE] - omega);
	//	dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE] 
	//					*tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] * (real(omega) -omega_star[ix][iL][iE])*Yp_R;
		dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE]
                                        *tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] * ((omega) -omega_star[ix][iL][iE])*Yp_R;
	}
	
	
	complex<double> sum;
	sum = simpintegral_3D(dwk_3D,grid->nx,grid->dr,grid->nL,grid->dL,grid->nE,grid->dE);
	return sum;
	
}


void help()
{
        cout<<"dwk++ is a code to calcualte delta wk and solver dispersion relation"<<endl;
        cout<<"the default input file is dwk.cfg"<<endl;
}

