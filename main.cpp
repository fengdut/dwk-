#include"tokamak.h"
#include"math.h"

#include<iostream>
#include<fstream>
#include <libconfig.h++>
#include"readinput.h"
#include"simpintegral.h"
#include"AllocArray.h"
#include"vector.h"




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

 	cout.precision (std::numeric_limits<long double>::digits10 + 1);
	cout<<"a: \t"<<tok.a<<"\tR0: \t"<<tok.R0<<"\teps: \t"<<tok.eps<<endl;
	cout<<"r0, rd: \t"<<slowing.r0<<"\t"<<slowing.rd<<endl;
	cout<<"L0, Ld: \t"<<slowing.L0<<"\t"<<slowing.Ld<<endl;
	cout<<"E0, Ed, Ec: \t"<<slowing.E0<<"\t"<<slowing.Ed<<"\t"<<slowing.Ec<<endl;	

	cout<<"nx, nL, nE, ntheta: \t"<<grid.nx<<"\t"<<grid.nL<<"\t"<<grid.nE<<"\t"<<grid.ntheta<<endl;
	cout<<"ra, rb: \t"<<grid.ra<<"\t"<<grid.rb<<endl;
	cout<<"La, Lb: \t"<<grid.La<<"\t"<<grid.Lb<<endl;
	cout<<"Ea, Eb: \t"<<grid.Ea<<"\t"<<grid.Eb<<endl;
	
//	cout<<"pi: \t"<<M_PI<<endl;
	
//to begin to calculate non-omega parts
	double *xarray,*Larray,*Earray;
	Alloc1D(xarray,grid.nx);
	Alloc1D(Larray,grid.nL);
	Alloc1D(Earray,grid.nE);
	linspace(grid.ra,grid.rb,grid.nx,xarray);	
	
	double *q_1D;
	Alloc1D(q_1D,grid.nx);
	qprofile(grid.nx,xarray,q_1D);
	//ofstream fout("test.dat");
	//printv(fout,grid.nx,q_1D);	
	v_max_min(grid.nx,q_1D);
	
	double ***F_3D;
	Alloc3D(F_3D,grid.nx,grid.nL,grid.nE);
		
	
	Free3D(F_3D);	
	Free1D(q_1D);
	Free1D(xarray);
	Free1D(Larray);
	Free1D(Earray);
	

}


