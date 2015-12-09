#ifndef VECTOR_H
#define VECTOR_H
#include<iostream>
#include<complex>
void linspace(const double xa,const double xb,const int nx,double *xarray);
void printv( std::ostream &fout,const int nx,const double *xarray);
void v_max_min(const int nx,double *xarray);
void max_min_2D(const int nx, const int ny, double ** F);
void max_min_2D(const int nx, const int ny, std::complex<double> ** F);
void max_min_3D(const int nx, const int ny, const int nz, double *** F);
#endif

