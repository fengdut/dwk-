#include <iostream>
#include<cmath>
#include<libconfig.h>
#include <libconfig.h++>
#include"tokamak.h"
#include"mode.h"
#include<complex>
#include<assert.h>
#include"AllocArray.h"
#include"readinput.h"
//#define EXIT_FAILURE 100
#include"outlog.h"

int read_tokamak(char* filename,Tokamak *ptok,Grid *pgrid,Slowing *pslowing,Mode *mode,Dwkopt *pdwkopt)
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
    		std::cerr << "Check the input file:" <<filename<<std::endl;
    		return(EXIT_FAILURE);
 	}
  	catch(const ParseException &pex)
 	{
    		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              	<< " - " << pex.getError() << std::endl;
    		std::cerr << "Check the input file:" <<filename<<std::endl;
    		return(EXIT_FAILURE);
 	}

/////////////////get tokamak parameters
	const Setting & root =cfg.getRoot();
	try
	{
		const Setting &tok =root["tokamak"];
		tok.lookupValue("a",ptok->a);
		tok.lookupValue("R0",ptok->R0);
		tok.lookupValue("Bt",ptok->Bt);
		tok.lookupValue("n0",ptok->n0);
		tok.lookupValue("mi",ptok->mi);
		tok.lookupValue("E_i0",ptok->E_i0);
		tok.lookupValue("m_ep",ptok->m_ep);
		
		Setting& set=tok["qc"];
		int lqc=set.getLength();
		int rlqc=0;
		if(lqc<nqc)
		{
			rlqc=lqc;
			for(int i=lqc;i<nqc;i++)
				ptok->qc[i]=0;
		}
		else
			rlqc=nqc;
		for(int i=0;i<rlqc;i++)
		{
			ptok->qc[i] = set[i];
		}
		
		tok.lookupValue("q_s",ptok->q_s);
	}
	catch(const SettingNotFoundException &nfex)
  	{
		cout<<"No tokamak parameters in filename \t"<<filename<<endl;
    		std::cerr << "Check the input file:" <<filename<<std::endl;
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
		assert(nx%3==1);
		assert(nL%3==1);
		assert(nE%3==1);
		assert(ntheta%3==1);
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
		int sigma;
		slowing.lookupValue("sigma",sigma);
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
		const Setting &modeset =root["mode"];
		int n,m,pa,pb;
		double r_s,delta_r;
		double omega_0,omega_1,omega_i;
		int omega_n;
		modeset.lookupValue("n",n);
		modeset.lookupValue("m",m);
		modeset.lookupValue("pa",pa);
		modeset.lookupValue("pb",pb);
		modeset.lookupValue("delta_r",delta_r);
		modeset.lookupValue("omega_0",omega_0);
		modeset.lookupValue("omega_1",omega_1);
		modeset.lookupValue("omega_i",omega_i);
		modeset.lookupValue("omega_n",omega_n);
		modeset.lookupValue("omega_err",mode->omega_err);
		modeset.lookupValue("max_iter",mode->max_iter);
		modeset.lookupValue("max_iterg",mode->max_iterg);
		modeset.lookupValue("dw_f",mode->dw_f);
		modeset.lookupValue("zero_rhod",mode->zero_rhod);
		modeset.lookupValue("xi_0",mode->xi_0);


		mode->n=n;
		mode->m=m;
		mode->pa =pa;
		mode->pb =pb;
		assert(pa<=pb);
		mode->delta_r=delta_r;
		mode->omega_0=omega_0;
		mode->omega_1=omega_1;
		mode->omega_i=omega_i;
		mode->omega_n=omega_n;
	
		assert(omega_n>1&&omega_n<10000);
		Alloc1D(mode->omega_array,omega_n);
		
		double domega=(omega_1-omega_0)/(omega_n-1);
		std::complex<double> ti=1.0i;
		for(int iomega=0;iomega<omega_n;iomega++)
		{
			mode->omega_array[iomega] =	omega_0 +domega*iomega +ti*omega_i; 
		}
	
	}
	 catch(const SettingNotFoundException &nfex)
        {
                cout<<"No mode parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
        }

	try
	{
		const Setting &optset =root["dwkopt"];
		optset.lookupValue("omega_star_off",pdwkopt->omega_star_off);
		optset.lookupValue("omega_off",pdwkopt->omega_off);

		
	}	
	catch(const SettingNotFoundException &nfex)
        {
                cerr<<"No mode parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
        }
	
	return 0;
}
