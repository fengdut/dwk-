#include"netcdf.h"
#include"stdio.h"
#include<iostream>
#include"stdlib.h"
#include"tokamak.h"
using namespace std;

#define ERRCODE -100
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int g_nr,g_nL,g_nE;


int open_netcdf(Grid *pgrid, char *filename)
{
	int ncid;
	//int err = nc_create(filename,NC_CLOBBER,&ncid);
	int err = nc_create(filename,NC_CLOBBER,&ncid);
	if(err)
	{
		cerr<<"open netcdf error"<<endl;
		exit(-1);
	}
	
//	cout<<"define dim 00"<<endl;
	err=nc_def_dim(ncid,"nr",pgrid->nx,&g_nr);
	if(err)
	 	ERR(err);
	
//	cout<<"define dim 01"<<endl;
	err=nc_def_dim(ncid,"nL",pgrid->nL,&g_nL);
	if(err)
         ERR(err);

//	cout<<"define dim 02"<<endl;
	err=nc_def_dim(ncid,"nE",pgrid->nE,&g_nE);
	if(err)
         ERR(err);

	int rid,Lid,Eid;
//	cout<<"define dim 03"<<endl;
	err=nc_def_var(ncid,"r",NC_DOUBLE,1,&g_nr,&rid);
	if(err)
         ERR(err);
	err=nc_def_var(ncid,"L",NC_DOUBLE,1,&g_nL,&Lid);
	if(err)
         ERR(err);
	err=nc_def_var(ncid,"E",NC_DOUBLE,1,&g_nE,&Eid);
	if(err)
         ERR(err);
	nc_enddef(ncid);

//	cout<<"define dim 04"<<endl;
	err=nc_put_var_double(ncid,rid,pgrid->xarray);	
	 if(err)
         ERR(err);
	err=nc_put_var_double(ncid,Lid,pgrid->Larray);	
	 if(err)
         ERR(err);
	err=nc_put_var_double(ncid,Eid,pgrid->Earray);	
	 if(err)
         ERR(err);
	return ncid;
			
}
void write_data_3D(double *** pdata, char * dataname,int fileid)
{
	static int init=0;
	static int dx,dy,dz;
	int err;

//	cout<<"define dim 0"<<endl;
	err =nc_redef(fileid);
	int dmins[3]={g_nr,g_nL,g_nE};
	int dataid;
//	cout<<"define dim 4"<<endl;
	nc_def_var(fileid,dataname,NC_DOUBLE,3,dmins,&dataid);
	
//	cout<<"define dim 5"<<endl;
	err=nc_enddef(fileid);

//	cout<<"define dim 6"<<endl;
	err=	nc_put_var_double(fileid,dataid,pdata[0][0]);
	if(err)
	ERR(err);
		
}
void write_data_2D(double **pdata, char *dataname,int fileid)
{
	
}
int close_netcdf(int fileid)
{
	nc_close(fileid);
}




