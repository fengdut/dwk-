#include"tokamak.h"
#include"math.h"

#include<iostream>
#include<fstream>
#include<ctime>
#include <libconfig.h++>
#include"readinput.h"
#include"simpintegral.h"
#include"AllocArray.h"
#include"vector.h"
#include"dist.h"




	using namespace std;
void help()
{
	cout<<"test delta wk"<<endl;
}

int main()
{	
	help();
	Tokamak tok;
	Grid 	grid;
	Slowing slowing;

	char filename[]="dwk.cfg";
//read input parameters
	read_tokamak(filename,&tok,&grid,&slowing);	

 	cout.precision (12);
	cout<<"a: \t"<<tok.a<<"\tR0: \t"<<tok.R0<<"\teps: \t"<<tok.eps<<endl;
	cout<<"r0, rd: \t"<<slowing.r0<<"\t"<<slowing.rd<<endl;
	cout<<"L0, Ld: \t"<<slowing.L0<<"\t"<<slowing.Ld<<endl;
	cout<<"E0, Ed, Ec: \t"<<slowing.E0<<"\t"<<slowing.Ed<<"\t"<<slowing.Ec<<endl;	

	cout<<"nx, nL, nE, ntheta: \t"<<grid.nx<<"\t"<<grid.nL<<"\t"<<grid.nE<<"\t"<<grid.ntheta<<endl;
	cout<<"ra, rb: \t"<<grid.ra<<"\t"<<grid.rb<<endl;
	cout<<"La, Lb: \t"<<grid.La<<"\t"<<grid.Lb<<endl;
	cout<<"Ea, Eb: \t"<<grid.Ea<<"\t"<<grid.Eb<<endl;
	
//	cout<<"pi: \t"<<M_PI<<endl;
	clock_t c_start=clock();
	
//to begin to calculate non-omega parts
	double *xarray,*Larray,*Earray, *tarray;
	Alloc1D(xarray,grid.nx);
	Alloc1D(Larray,grid.nL);
	Alloc1D(Earray,grid.nE);
	Alloc1D(tarray,grid.ntheta);
	linspace(grid.ra,grid.rb,grid.nx,xarray);	
	linspace(grid.La,grid.Lb,grid.nL,Larray);
	linspace(grid.Ea,grid.Eb,grid.nE,Earray);
	linspace(0,2*M_PI,grid.ntheta,tarray);
	grid.xarray=xarray;
	grid.Larray=Larray;
	grid.Earray=Earray;
	grid.tarray=tarray;
	
	double *q_1D;
	Alloc1D(q_1D,grid.nx);
	qprofile(grid.nx,xarray,q_1D);
	//ofstream fout("test.dat");
	//printv(fout,grid.nx,q_1D);	
	v_max_min(grid.nx,q_1D);
	
	double ***F_3D;
	Alloc3D(F_3D,grid.nx,grid.nL,grid.nE);
		
	F0_3D_test(&slowing,&grid,F_3D);	
	Free3D(F_3D);	
	Free1D(q_1D);
	Free1D(xarray);
	Free1D(Larray);
	Free1D(Earray);
	
	cout<<"time:"<<float(clock() - c_start)/CLOCKS_PER_SEC<<endl;
}


