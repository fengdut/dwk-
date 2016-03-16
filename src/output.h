#include"netcdf.h"
#include<iostream>
#include"stdlib.h"
#include<fstream>
using namespace std;



int open_netcdf(Grid *pgrid, char *filename);
void write_data_3D(double *** pdata, char * dataname,int fileid);
int close_netcdf(int fileid);


