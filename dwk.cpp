#include"tokamak.h"
#include<complex>
#include<cmath>
#include<iostream>
#include<fstream>
#include <algorithm>
#include"AllocArray.h"
#include"simpintegral.h"
#include"mode.h"
#include"particle_fre.h"
#include"vector.h"

using namespace std;


void help()
{
	cout<<endl;
        cout<<"dwk++ is a code to calcualte delta wk for energetic particle"<<endl;
	cout<<"in tokamak, and solver dispersion relation."<<endl;
	cout<<endl;
        cout<<"-i <inputfile>, the default input file is dwk.cfg."<<endl;
	cout<<"-o <outputfile>, the default output file is dwk_omega_dwk.out."<<endl;
	cout<<"-s scan dwk(omega). Or only find omega_0 and beta_c."<<endl;
	cout<<endl;
	
}


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

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
//		complex<double> Yp_R =Yp2_2D[ix][iL] /(0.5 +p *omega_b_3D[ix][iL][iE] - omega);
		dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE] 
						*tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] * (real(omega) -omega_star[ix][iL][iE])*Yp_R;
//		dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE]
                                  //     *tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] * (0.0*(omega) -omega_star[ix][iL][iE])*Yp_R; //omega_star only
	                //dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*(grid->Earray[iE])
                         //               *tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] * (real(omega) -0.0*omega_star[ix][iL][iE])*Yp_R; //omega only
	}
	
	
	complex<double> sum;
	sum = simpintegral_3D(dwk_3D,grid->nx,grid->dr,grid->nL,grid->dL,grid->nE,grid->dE);
	return sum;
	
}

//dwk(omega)
complex<double> dwk_omega(Grid *const grid,Mode *const mode,complex<double> omega,
        double *** const omega_phi_3D, double ***const omega_b_3D, double *** tau_b_3D,
        complex<double> ** const Yps_2D, double * const J_q_1D, double *** const F_E_3D,
        double *** const omega_star,
        complex<double> **const G_2D,double **const Chi_2D,double *** const b_lambda_3D, double *** const lambda_b_3D,
        double ***const Theta_3D)
{

        complex<double> ***dwk_3D;
	complex<double>dwk=0;
        Alloc3D(dwk_3D,grid->nx,grid->nL,grid->nE);
        for(int p=mode->pa;p<=mode->pb;p++)
        {
                Yps(grid,G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,p,Yps_2D);
		//cout<<"Yps_2D ";
		//max_min_2D( grid->nx,  grid->nL, Yps_2D);
                for(int ix=0;ix<grid->nx;ix++)
                for(int iL=0;iL<grid->nL;iL++)
                for(int iE=0;iE<grid->nE;iE++)
                {
                        complex<double> Yp_R =abs(Yps_2D[ix][iL])*abs(Yps_2D[ix][iL])
                        		/(mode->n*omega_phi_3D[ix][iL][iE] +p *omega_b_3D[ix][iL][iE] - omega);
                        dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE]
                        		*tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE]
                        		*(real(omega) -omega_star[ix][iL][iE])*Yp_R;
                }
                dwk += simpintegral_3D(dwk_3D,grid->nx,grid->dr,grid->nL,grid->dL,grid->nE,grid->dE);
        }
        Free3D(dwk_3D);
	return dwk;
}

//dwk(omega)
void dwk_omega_array(Grid *const grid,Mode *const mode, 
        double *** const omega_phi_3D, double ***const omega_b_3D, double *** tau_b_3D,
        complex<double> ** const Yps_2D, double * const J_q_1D, double *** const F_E_3D, 
	double *** const omega_star, 
	complex<double> **const G_2D,double **const Chi_2D,double *** const b_lambda_3D, double *** const lambda_b_3D,
	double ***const Theta_3D,
	complex<double> * dwk_array,char const *filename)
{
	
	cout<<"------omega vs. dwk-------"<<endl;
	complex<double> ***dwk_3D;
	Alloc3D(dwk_3D,grid->nx,grid->nL,grid->nE);
	cout.precision (12);
	
	ofstream fout(filename);	
	fout.precision (12);
	for(int iw=0;iw<mode->omega_n;iw++)
	{
		complex<double> omega =mode->omega_array[iw];
		for(int p=mode->pa;p<=mode->pb;p++)
		{

			Yps(grid,G_2D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,p,Yps_2D);
	               //cout<<"Yps_2D ";
                //	max_min_2D( grid->nx,  grid->nL, Yps_2D);
	         //      cout<<"G_2D ";
                //	max_min_2D( grid->nx,  grid->nL, G_2D);
			for(int ix=0;ix<grid->nx;ix++)
        		for(int iL=0;iL<grid->nL;iL++)
        		for(int iE=0;iE<grid->nE;iE++)
        		{
                		complex<double> Yp_R =abs(Yps_2D[ix][iL])*abs(Yps_2D[ix][iL]) 
						/(mode->n*omega_phi_3D[ix][iL][iE] +p *omega_b_3D[ix][iL][iE] - omega);
                		dwk_3D[ix][iL][iE] = J_q_1D[ix] *grid->Earray[iE]*grid->Earray[iE]*grid->Earray[iE]
                                                     *tau_b_3D[ix][iL][iE] *F_E_3D[ix][iL][iE] 
						     *(real(omega) -omega_star[ix][iL][iE])*Yp_R;
        		}
			dwk_array[iw] += simpintegral_3D(dwk_3D,grid->nx,grid->dr,grid->nL,grid->dL,grid->nE,grid->dE);
		}
		cout<<real(omega)<<'\t'<<imag(omega)<<'\t'<<real(dwk_array[iw])<<'\t'<<imag(dwk_array[iw])<<endl;
		fout<<real(omega)<<'\t'<<imag(omega)<<'\t'<<real(dwk_array[iw])<<'\t'<<imag(dwk_array[iw])<<endl;
	//	cout<<real(omega)<<'\t'<<endl;
	}
		
	cout<<"-----end omega vs. dwk--------"<<endl;
	Free3D(dwk_3D);
}

//find the real(omega_0) where real(dwk(omega_0))==0;
complex<double> find_dwk_omega0(Grid *const grid,Mode *const mode,Tokamak *tok,
        double *** const omega_phi_3D, double ***const omega_b_3D, double *** tau_b_3D,
        complex<double> ** const Yps_2D, double * const J_q_1D, double *** const F_E_3D,
        double *** const omega_star,
        complex<double> **const G_2D,double **const Chi_2D,double *** const b_lambda_3D, double *** const lambda_b_3D,
        double ***const Theta_3D,complex<double> *dwk_0)
{
        cout<<"------------------- find omega_0 ----------------------"<<endl;
        complex<double> ***dwk_3D;
        Alloc3D(dwk_3D,grid->nx,grid->nL,grid->nE);
        cout.precision (12);
        complex<double> omega_a =mode->omega_array[0];
	complex<double> omega_b =mode->omega_array[mode->omega_n-1];
	complex<double> dwk_a,dwk_b,dwk_test,omega_0;

	dwk_a= 	dwk_omega(grid,mode,/* */ omega_a /**/, omega_phi_3D,omega_b_3D, tau_b_3D, Yps_2D, J_q_1D, F_E_3D, omega_star,
        	G_2D, Chi_2D, b_lambda_3D, lambda_b_3D, Theta_3D);
	dwk_b= 	dwk_omega(grid,mode,/* */ omega_b /**/, omega_phi_3D,omega_b_3D, tau_b_3D, Yps_2D, J_q_1D, F_E_3D, omega_star,
        	G_2D, Chi_2D, b_lambda_3D, lambda_b_3D, Theta_3D);
        cout<<"omega_a: "<<real(omega_a)<<'\t'<<imag(omega_a)<<'\t'<<real(dwk_a)<<'\t'<<imag(dwk_a)<<endl;
        cout<<"omega_b: "<<real(omega_b)<<'\t'<<imag(omega_b)<<'\t'<<real(dwk_b)<<'\t'<<imag(dwk_b)<<endl;
	
	double gamma=0;
	double Cb=1;
	if(tok->beta_h==0)
	{
		gamma=0.0;
		Cb=1.0;
	}
	else
	{
		gamma=imag(omega_a);
		Cb =tok->C*tok->beta_h;
	}
	cout<<"******************"<<endl;
	cout<<"gamma: "<<gamma<<endl;
	cout<<"******************"<<endl;
	
	if((Cb*real(dwk_a)+gamma)*(Cb*real(dwk_b)+gamma)<=0)
	{
		int in=0;
		while((Cb*real(dwk_a)+gamma)*(Cb*real(dwk_b)+gamma)<=0)
		{	
			complex<double> omega_test = 0.5*(omega_a+omega_b);
			dwk_test =dwk_omega(grid,mode,/* */ omega_test /**/,
				 omega_phi_3D,omega_b_3D, tau_b_3D, Yps_2D, J_q_1D, F_E_3D, omega_star,
                	         G_2D, Chi_2D, b_lambda_3D, lambda_b_3D, Theta_3D); 
			if(real(omega_b)-real(omega_a)<mode->omega_err)	
			{
				omega_0=omega_test;
				cout<<"omega_0: "<<real(omega_test)<<'\t'<<imag(omega_test)
				    <<"\t dwk_0: "<<real(dwk_test)<<'\t'<<imag(dwk_test)<<endl;
				break;
			}
			else
			{
				if((Cb*real(dwk_a)+gamma)*(Cb*real(dwk_test)+gamma)<=0)
				{
					omega_b=omega_test;
					dwk_b  =dwk_test;
				}
				else
				{
					omega_a=omega_test;
					dwk_a = dwk_test;
				}
				cout<<in<<"\t omega_a"<<omega_a <<"\t omega_b" <<omega_b<<"\terr:"<<real(omega_b-omega_a)<<endl; 
			}
			in++;
			if(in>mode->max_iter)
			{
				cout<<"max iter reach"<<endl;
				break;
			}
			
		}
	}
	else
	{
		cout<<"*************change the omega range*****************"<<endl;
	}
	
        cout<<"-----end find omega_0--------"<<endl;
        Free3D(dwk_3D);
	*dwk_0 = dwk_test;
	return omega_0;
}





