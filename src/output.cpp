#include"netcdf.h"
#include<iostream>
#include"stdlib.h"
using namespace std;

#define ERRCODE -100
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
int open_netcdf(char *filename)
{
	int ncid;
	int err = nc_create(filename,NC_CLOBBER,&ncid);
	if(err)
	{
		cerr<<"open netcdf error"<<endl;
		exit(-1);
	}
	return ncid;
			
}
void write_data_3D(int nx,int ny, int nz, double *** pdata, char * dataname,int fileid)
{
	int err;
	int dx,dy,dz;
	err=	nc_def_dim(fileid,"nx",nx,&dx);
     	if(err)
        	ERR(err);

	nc_def_dim(fileid,"ny",ny,&dy);
	nc_def_dim(fileid,"nz",nz,&dz);

	int dmins[3]={dx,dy,dz};
	int dataid;
	nc_def_var(fileid,dataname,NC_DOUBLE,3,dmins,&dataid);

	err=	nc_put_var_double(fileid,dataid,pdata[0][0]);
	if(err)
	ERR(err);
		
}

int close_netcdf(int fileid)
{
	nc_close(fileid);
}




