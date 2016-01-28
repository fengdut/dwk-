#include <iostream>
#include<libconfig.h>
#include <libconfig.h++>
#include"tokamak.h"

#define EXIT_FAILURE 100

int read_tokamak(char* filename,Tokamak *ptok,Grid *pgrid,Slowing *pslowing)
{
	using namespace std;
	using namespace libconfig;

	Config cfg;
///////////open config file
	try
	{
		cfg.readFile(filename);
	}
	catch(const FileIOException &fioex)
  	{
    		std::cerr << "I/O error while reading file." << std::endl;
    		return(EXIT_FAILURE);
 	}
  	catch(const ParseException &pex)
 	{
    		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              	<< " - " << pex.getError() << std::endl;
    		return(EXIT_FAILURE);
 	}

/////////////////get tokamak parameters
	const Setting & root =cfg.getRoot();
	try
	{
		const Setting &tok =root["tokamak"];
		double a;
		double R0;
		tok.lookupValue("a",a);
		tok.lookupValue("R0",R0);
		ptok->a=a;
		ptok->R0=R0;
		ptok->eps=a/R0;
	}
	catch(const SettingNotFoundException &nfex)
  	{
		cout<<"No tokamak parameters in filename \t"<<filename<<endl;
		return(EXIT_FAILURE);
  	}
////////////////get mesh config
	 try
        {
                const Setting &grid =root["grid"];
		int nx,nL,nE,ntheta;
                grid.lookupValue("nx",nx);
                grid.lookupValue("nL",nL);
                grid.lookupValue("nE",nE);
                grid.lookupValue("ntheta",ntheta);
		pgrid->nx =nx;
		pgrid->nL =nL;
		pgrid->nE =nE;
		pgrid->ntheta=ntheta;
        }
        catch(const SettingNotFoundException &nfex)
        {
                cout<<"No grid parameters in filename \t"<<filename<<endl;
		return(EXIT_FAILURE);
        }
/////////////////get slowing down distribution parameters
	try
	{
		const Setting &slowing =root["slowing"];	
		double rd,r0,L0,Ld,E0,Ed,Ec;
		slowing.lookupValue("r0",r0);	
		slowing.lookupValue("rd",rd);	
		slowing.lookupValue("L0",L0);	
		slowing.lookupValue("Ld",Ld);	
		slowing.lookupValue("E0",E0);	
		slowing.lookupValue("Ed",Ed);	
		slowing.lookupValue("Ec",Ec);	
		pslowing->r0 =r0;
		pslowing->rd =rd;
		pslowing->L0 =L0;
		pslowing->Ld =Ld;
		pslowing->E0 =E0;
		pslowing->Ed =Ed;
		pslowing->Ec =Ec;
	}	
	catch(const SettingNotFoundException &nfex)
	{
		cout<<"No slowing down parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
	}
	
	pgrid->ra=1e-6;
	pgrid->rb=1;
	pgrid->dr= (pgrid->rb- pgrid->ra)/(pgrid->nx-1);

	pgrid->La=1e-6;
	pgrid->Lb=0.1;
	pgrid->dL= (pgrid->Lb - pgrid->La)/(pgrid->nL-1);

	pgrid->Ea=0.1;
	pgrid->Eb=pslowing->E0;	
	pgrid->dE= (pgrid->Eb -pgrid->Ea)/(pgrid->nE-1);	
	
	return 0;
}
