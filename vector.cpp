#ifndef VECTOR_H
#define VECTOR_H

#include<iostream>
#include<complex>
	using namespace std;
void linspace(const double xa,const double xb,const int nx,double *xarray)
{
	double dx=(xb-xa)/(nx-1);

	for(int i=0;i<nx;i++)
	{
		xarray[i] =xa +i*dx;		
	}
}

void printv( std::ostream &fout,const int nx,const double *xarray)
{
using namespace std;
	for(int i=0;i<nx;i++)
	{
		fout<<i<<"\t"<<xarray[i]<<endl;
	}
}

void v_max_min(const int nx,double *xarray)
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

void max_min_2D(const int nx, const int ny, double ** F)
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

void max_min_3D(const int nx, const int ny, const int nz, double *** F)
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

void max_min_2D(const int nx, const int ny, complex<double> ** F)
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

#endif

