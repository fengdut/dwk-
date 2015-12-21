#ifndef VECTOR_H
#define VECTOR_H
#include<iostream>
#include<complex>
void linspace(double const xa,const double xb,int const nx,double *xarray);
void printv( std::ostream &fout,int const nx,double *const xarray);
void max_min_1D(int const nx,double *const xarray);
void max_min_2D(int const nx, const int ny, double **const  F);
void max_min_2D(int const nx, const int ny, std::complex<double> **const  F);
void max_min_3D(int const nx, const int ny, const int nz, double ***const F);
void max_min_3D(int const nx, int const  ny, int const nz,std::complex<double>*** const  F);
#endif

