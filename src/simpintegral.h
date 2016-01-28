#ifndef SIMPINTEGRAL_H
#define SIMPINTEGRAL_H

#include"AllocArray.h"
#include<iostream>
#include<assert.h>

//Simpson's 3/8 rule
//ny must be equiv to 3*n+1
//double simpintegral(double *y,int ny,double dx);

template <class datatype>
datatype simpintegral(datatype * const y, int const  ny, double const dx)
{
	assert(ny%3==1);
        datatype sum=0;
        for(int i=0;i<ny;i++)
        {
                sum +=3.0*y[i];
        }
        for(int i=3;i<ny-1;i=i+3)
        {
                sum -= y[i];
        }
        sum =sum -2.0*(y[0]+y[ny-1]);
        sum =dx*3.0*sum*0.125;
        return sum;
}

template <class datatype>
datatype simpintegral_o(datatype * const y, int const  ny, double const dx)
{
        datatype sum=0;
	float n0=0;
        for(int i=0;i<ny;i++)
        {
		float d=floor(i/3)-n0;
                sum +=(3.0-d)*y[i];
		n0=n0+d;
        }
        sum =sum -2.0*(y[0]+y[ny-1]);
        sum =dx*3.0*sum*0.125;
        return sum;
}

template <class datatype>
datatype simpintegral_2D(datatype** const F_2D, int const nx, double const dx, int const ny,double const dy)
{
        datatype sum =0;
        datatype *F_1D;
        Alloc1D(F_1D,nx);
	using namespace std;

        for(int ix=0;ix<nx;ix++)
	{
                F_1D[ix] = simpintegral(F_2D[ix],ny,dy);
	//	cout<<"s_1: "<<F_1D[ix]<<endl;
	}
        sum = simpintegral(F_1D,nx,dx);

        Free1D(F_1D);
        return sum;
}



template <class datatype>
datatype simpintegral_3D(datatype*** const F_3D, int const nx, double const dx, int const ny,double const dy, int const nz, double const dz)
{
	datatype sum =0;
	datatype ** F_2D,*F_1D;
	Alloc2D(F_2D,nx,ny);
	Alloc1D(F_1D,nx);
	for(int ix=0;ix<nx;ix++)
	for(int iy=0;iy<ny;iy++)
	{
		F_2D[ix][iy] = simpintegral(F_3D[ix][iy],nz,dz);
	}
	
	for(int ix=0;ix<nx;ix++)
		F_1D[ix] = simpintegral(F_2D[ix],ny,dy);
	
	sum = simpintegral(F_1D,nx,dx);
	
	Free2D(F_2D);
	Free1D(F_1D);
	return sum;	
}

#endif

