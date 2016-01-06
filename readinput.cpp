#include <iostream>
#include<cmath>
#include<libconfig.h>
#include <libconfig.h++>
#include"tokamak.h"
#include"mode.h"


#define EXIT_FAILURE 100

int read_tokamak(char* filename,Tokamak *ptok,Grid *pgrid,Slowing *pslowing,Mode *mode)
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
		double C;
		tok.lookupValue("a",a);
		tok.lookupValue("R0",R0);
		tok.lookupValue("C",C);

		ptok->a=a;
		ptok->R0=R0;
		ptok->eps=a/R0;
		ptok->C=C;
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
		double rho_h,rho_d;
		int sigma;
		slowing.lookupValue("rho_h",rho_h);
		slowing.lookupValue("rho_d",rho_d);
		slowing.lookupValue("sigma",sigma);
		pslowing->rho_h = rho_h;
		pslowing->rho_d = rho_d;
		pslowing->sigma = sigma;
	}	
	catch(const SettingNotFoundException &nfex)
	{
		cout<<"No slowing down parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
	}
//////////// get mode parameters  
	try
	{
		const Setting &slowing =root["mode"];
		int n,m,pa,pb;
		double r_s,delta_r;
		slowing.lookupValue("n",n);
		slowing.lookupValue("m",m);
		slowing.lookupValue("pa",pa);
		slowing.lookupValue("pb",pb);
		slowing.lookupValue("r_s",r_s);
		slowing.lookupValue("delta_r",delta_r);
		mode->n=n;
		mode->m=m;
		mode->pa =pa;
		mode->pb =pb;
		mode->r_s=r_s;
		mode->delta_r=delta_r;
	}
	catch(const SettingNotFoundException &nfex)
        {
                cout<<"No mode parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
        }



	return 0;
}
