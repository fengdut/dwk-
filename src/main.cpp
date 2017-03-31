#include<fstream>
#include<iostream>
#include<math.h>
#include<ctime>
#include<complex>

#include"vector.h"
#include"AllocArray.h"
#include"tokamak.h"
#include"readinput.h"
#include"dist.h"
#include"mode.h"
#include"particle_fre.h"
#include"dwk.h"
#include"output.h"
#include"outlog.h"
#include"sys/time.h"

outlog log_cout("dwk.log");


using namespace std;
int main(int arg,char * argx[])
{	
	struct timeval start, stop, diff;
	gettimeofday(&start,0);
		
	if(cmdOptionExists(argx, argx+arg, "-h"))//print help information only
	{
		help();
		return 0;
	}
	char default_ifilename[]="dwk.cfg";
	char default_ofilename[]="omega_dwk.out";
	char *ifilename =getCmdOption(argx,argx+arg,"-i");
	char *ofilename =getCmdOption(argx,argx+arg,"-o");

	if(ifilename==0)
		ifilename=default_ifilename;
        if(ofilename==0)
                ofilename=default_ofilename;
	
	Tokamak tok;
	Grid 	grid;
	Slowing slowing;
	Mode mode;
	Dwkopt dwkopt;
//read input parameters
	read_tokamak(ifilename,&tok,&grid,&slowing,&mode,&dwkopt);	
	CGrid pgrid(&grid,&slowing);
 	cout.precision (8);
	pgrid.showgrid();	
//alloc memory
	double *q_1D,*J_q_1D;
	Alloc1D(q_1D,grid.nx);
	Alloc1D(J_q_1D,grid.nx);
	double **Chi_2D, **K_2D,**kappa_2D;	
	complex<double> ***G_3D;
	Alloc3D(G_3D,grid.nx,grid.nE,grid.ntheta);

	Alloc2D(Chi_2D,grid.nx,grid.nL);
	Alloc2D(K_2D,grid.nx,grid.nL);
	Alloc2D(kappa_2D,grid.nx,grid.nL);
	double *** omega_b_3D,***omega_phi_3D,*** tau_b_3D;
	Alloc3D(omega_b_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_phi_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(tau_b_3D,grid.nx,grid.nL,grid.nE);
	double ***F_3D,***FE_3D,***FR_3D, ***omega_star_3D;
	Alloc3D(F_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FE_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(FR_3D,grid.nx,grid.nL,grid.nE);
	Alloc3D(omega_star_3D,grid.nx,grid.nL,grid.nE);
	double *** lambda_b_3D, ***b_lambda_3D, ***Theta_3D;
	Alloc3D(lambda_b_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(b_lambda_3D,grid.nx,grid.nL,grid.ntheta);
	Alloc3D(Theta_3D,grid.nx,grid.nL,grid.ntheta);
	complex<double> ***dwk_3D;
	Alloc3D(dwk_3D,grid.nx,grid.nL,grid.nE);
//end alloc memory

	       gettimeofday(&stop,0);
        cout<<"1111 dwk++ Elapsed time:\t"<<1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec-start.tv_usec<<endl;

//begin calculate non-omega parts------------------------
	qprofile(&grid,&tok,q_1D);
	find_rs(&grid,q_1D, &tok);
	calculate_normalization(&tok, &slowing,&mode);
	showtokamak(&tok,&slowing);
	double Cbeta=0;
	F0_3D(&slowing,&grid,&tok,slowing.rho_h,q_1D,mode.n,F_3D,FE_3D,FR_3D,omega_star_3D,&Cbeta);	
       gettimeofday(&stop,0);
        cout<<"1.2dwk++ Elapsed time:\t"<<1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec-start.tv_usec<<endl;
	Lambda_b_L_3D(&grid,&tok,lambda_b_3D,b_lambda_3D);
	Theta(b_lambda_3D,&grid,Theta_3D);
	G_R_theta(&grid,&tok,&slowing,&mode,q_1D,G_3D); 
	Chi(&grid,&tok,slowing.sigma,Chi_2D,kappa_2D,K_2D);	
	omega_b(&grid, &tok,kappa_2D, K_2D, q_1D,omega_b_3D);
	omega_phi(&grid,q_1D,omega_b_3D,omega_phi_3D);
	       gettimeofday(&stop,0);
        cout<<"2222 dwk++ Elapsed time:\t"<<1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec-start.tv_usec<<endl;
	tau_b(&grid,omega_b_3D,tau_b_3D);
	J_q(&grid,q_1D,J_q_1D);
//end calculate non-omega parts--------------------------	

	complex<double> *dwk_array;
	Alloc1D(dwk_array,mode.omega_n);
	
	if(cmdOptionExists(argx,argx+arg,"-s"))
        {
		cout<<"*********************************"<<endl;
		cout<<"scan dwk(omega)"<<endl;
                dwk_omega_array(&grid, &mode, omega_phi_3D, omega_b_3D, tau_b_3D,
                        J_q_1D, FE_3D, omega_star_3D,
                        G_3D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,
                        dwk_array,ofilename,&dwkopt);
		cout<<"outputfile: "<<ofilename<<endl;
		cout<<"end scan"<<endl;
		cout<<"*********************************"<<endl;
        }

	cout<<"*********************************"<<endl;
	if(cmdOptionExists(argx,argx+arg,"-y"))
        {
                cout<<"write Yps_3D to Yps.nc"<<endl;
                double *** rYps,***iYps;
                Alloc3D(rYps,grid.nx,grid.nL,grid.nE);
                Alloc3D(iYps,grid.nx,grid.nL,grid.nE);

                int fileid=0;
                char filename[]="Yps.nc";
                fileid =open_netcdf(&grid,filename);
                char dataname[10];

                int np=mode.pb-mode.pa+1;
                for(int i=0;i<np;i++)
                {
                for(int ix=0;ix<grid.nx;ix++)
                for(int iL=0;iL<grid.nL;iL++)
                for(int iE=0;iE<grid.nE;iE++)
                {
                        rYps[ix][iL][iE] = real(gYps_3D[i][ix][iL][iE]);
                        iYps[ix][iL][iE] = imag(gYps_3D[i][ix][iL][iE]);

                }
                sprintf(dataname,"rYps_%d",mode.pa+i);
                write_data_3D(rYps,dataname,fileid);
                sprintf(dataname,"iYps_%d",mode.pa+i);
                write_data_3D(iYps,dataname,fileid);
                }
                close_netcdf(fileid);
                Free3D(rYps);
                Free3D(iYps);

        }
       gettimeofday(&stop,0);
        cout<<"dwk++ Elapsed time:\t"<<1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec-start.tv_usec<<endl;

	if(!cmdOptionExists(argx,argx+arg,"-x"))
	{
	cout<<"The test run, assume imag(omega)=0 "<<endl;
	complex<double> omega_0,dwk_0;
	omega_0= find_dwk_omega0(&grid,&mode,&tok,omega_phi_3D,omega_b_3D,tau_b_3D,
			J_q_1D,FE_3D,omega_star_3D,
			G_3D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,&dwk_0,&dwkopt);
	tok.beta_h = real(omega_0)/(imag(dwk_0)*tok.C);
	cout<<"beta_h_0= "<<tok.beta_h*Cbeta<<endl;
	complex<double> err=0;
	complex<double> ti =1.0i;
	err =ti*omega_0 -tok.C*tok.beta_h*dwk_0;
	cout<<scientific;
	cout<<"error=i*omega_0 -C*beta_h *dwk(omega_0)= "<<err<<endl;
	cout<<fixed;
	cout<<"*********************************"<<endl;
	
	for(int gi=0;gi<mode.max_iterg;gi++)
	{
		cout<<"*********************************"<<endl;
		cout<<"find the solution"<<endl;
		omega_0= find_dwk_omega0(&grid,&mode,&tok,omega_phi_3D,omega_b_3D,tau_b_3D,
                        J_q_1D,FE_3D,omega_star_3D,
                        G_3D,Chi_2D,b_lambda_3D,lambda_b_3D,Theta_3D,&dwk_0,&dwkopt);
        	tok.beta_h = real(omega_0)/(imag(dwk_0)*tok.C);
        	cout<<"beta_h_"<<gi+1<<"= "<<tok.beta_h*Cbeta<<endl;
		err =ti*omega_0 -tok.C*tok.beta_h*dwk_0;
		cout<<scientific;
        	cout<<"error=i*omega_0 -C*beta_h *dwk(omega_0)= "<<err<<endl;
		cout<<fixed;
        	cout<<"*********************************"<<endl;
		if(abs(err)<mode.omega_err*tok.C*tok.beta_h)
			break;
	}
	
	cout<<"Omega_0= "<<real(omega_0)<<"+"<<imag(omega_0)<<"i\t"<<"dwk= "<<real(dwk_0)<<"+"<<imag(dwk_0)<<endl;
	cout<<"Beta_h= "<<tok.beta_h*Cbeta<<endl;
	cout <<"in omega_A unit"<<endl;
	cout <<"Omega_0A="<<omega_0*tok.omega_i0/tok.omega_A<<endl;
	cout <<"in kHz " <<endl;
	cout <<"Omega_0kHz="<<omega_0*tok.omega_i0/2.0/M_PI/1000.0 <<endl;
	cout<<"*************************************************************"<<endl;
	}
	if(cmdOptionExists(argx,argx+arg,"-o"))
	{
		cout <<"write omega_phi_3D to omega_phi.nc"<<endl;
		int fileid=0;
                char filename[]="omega_phi.nc";
                fileid =open_netcdf(&grid,filename);
                char dataname[10];
		sprintf(dataname,"omega_phi_3D_x_L_E");
		write_data_3D(omega_phi_3D,dataname,fileid);
	        close_netcdf(fileid);

	}
	Free1D(dwk_array);	
	Free1D(J_q_1D);		Free1D(q_1D);
        Free2D(Chi_2D); 	Free3D(G_3D);
        Free2D(kappa_2D);	Free2D(K_2D);
	Free3D(dwk_3D); 	Free3D(tau_b_3D);
	Free3D(omega_phi_3D);	Free3D(omega_b_3D);
	Free3D(Theta_3D);
	Free3D(b_lambda_3D);	Free3D(lambda_b_3D);
	Free3D(F_3D);		Free3D(FR_3D);	
	Free3D(FE_3D);		Free3D(omega_star_3D);
	gettimeofday(&stop,0);
	cout<<"dwk++ Elapsed time:\t"<<1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec-start.tv_usec<<endl;
	//cout<<"dwk++ Elapsed time:\t"<<float(clock() - c_start)/CLOCKS_PER_SEC<<" second"<<endl;
	fflush(stdout);
}


