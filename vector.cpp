#ifndef VECTOR_H
#define VECTOR_H

#include<iostream>
#include<complex>
	using namespace std;
void linspace(double const xa,const double xb, int const nx,double *xarray)
{
	double dx=(xb-xa)/(nx-1);

	for(int i=0;i<nx;i++)
	{
		xarray[i] =xa +i*dx;		
	}
}

void printv( std::ostream &fout, int const nx,double *const xarray)
{
using namespace std;
	for(int i=0;i<nx;i++)
	{
		fout<<i<<"\t"<<xarray[i]<<endl;
	}
}

void max_min_1D( int const nx,double *const xarray)
{
	double maxf,minf;
	maxf=minf=xarray[0];
        for(int i=1;i<nx;i++)
        {
		maxf =maxf>xarray[i] ? maxf:xarray[i];
		minf =minf<xarray[i] ? minf:xarray[i];
        }
	cout<<"max= "<<maxf<<"\tmin="<<minf<<endl;	
}

void max_min_2D( int const nx,  int const  ny, double **const F)
{
        double maxf,minf;
        maxf=minf=F[0][0];
        for(int ix=0;ix<nx;ix++)
        for(int iy=0;iy<ny;iy++)
        {
                maxf =maxf>F[ix][iy] ? maxf :F[ix][iy];
                minf =minf<F[ix][iy] ? minf :F[ix][iy];
        }
        cout<<"max= "<<maxf<<"\t min= "<<minf<<endl;


}

void max_min_3D( int const nx,  int const  ny, int const nz, double*** const F)
{
	double maxf,minf;
	maxf=minf=F[0][0][0];
	for(int ix=0;ix<nx;ix++)
	for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++)
	{
		maxf =maxf>F[ix][iy][iz] ? maxf :F[ix][iy][iz];
		minf =minf<F[ix][iy][iz] ? minf :F[ix][iy][iz];
	}
	cout<<"max= "<<maxf<<"\t min= "<<minf<<endl;

	
}

void max_min_2D( int const nx,  int const  ny, complex<double>** const  F)
{
        double maxfr,minfr;
        double maxfi,minfi;
        maxfr=minfr=real(F[0][0]);
        maxfi=minfi=imag(F[0][0]);
        for(int ix=0;ix<nx;ix++)
        for(int iy=0;iy<ny;iy++)
        {
                maxfr =maxfr>real(F[ix][iy]) ? maxfr :real(F[ix][iy]);
                minfr =minfr<real(F[ix][iy]) ? minfr :real(F[ix][iy]);
                maxfi =maxfi>imag(F[ix][iy]) ? maxfi :imag(F[ix][iy]);
                minfi =minfi<imag(F[ix][iy]) ? minfi :imag(F[ix][iy]);
        }
        cout<<"max(real)= "<<maxfr<<"\t max(imag)"<<maxfi<<"\t min(real)= "<<minfr<<"\t min(imag)"<<minfi<<endl;


}

void max_min_3D(int const nx, int const  ny, int const nz,complex<double>*** const  F)
{
        double maxfr,minfr;
        double maxfi,minfi;
        maxfr=minfr=real(F[0][0][0]);
        maxfi=minfi=imag(F[0][0][0]);
        for(int ix=0;ix<nx;ix++)
        for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++)
        {
                maxfr =maxfr>real(F[ix][iy][iz]) ? maxfr :real(F[ix][iy][iz]);
                minfr =minfr<real(F[ix][iy][iz]) ? minfr :real(F[ix][iy][iz]);
                maxfi =maxfi>imag(F[ix][iy][iz]) ? maxfi :imag(F[ix][iy][iz]);
                minfi =minfi<imag(F[ix][iy][iz]) ? minfi :imag(F[ix][iy][iz]);
        }
        cout<<"max(real)= "<<maxfr<<"\t max(imag)"<<maxfi<<"\t min(real)= "<<minfr<<"\t min(imag)"<<minfi<<endl;


}


#endif

