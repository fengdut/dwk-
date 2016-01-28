#ifndef VECTOR_H
#define VECTOR_H

#include<iostream>
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
	using namespace std;
	double maxf,minf;
	maxf=minf=xarray[0];
        for(int i=1;i<nx;i++)
        {
		maxf =maxf>xarray[i] ? maxf:xarray[i];
		minf =minf<xarray[i] ? minf:xarray[i];
        }
	cout<<"max, min: \t"<<maxf<<"\t"<<minf<<endl;	
}

#endif

