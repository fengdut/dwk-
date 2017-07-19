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
#define EXIT_FAILURE 100
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
		tok.lookupValue("beta_h",ptok->beta_h);
		tok.lookupValue("beta_hb",ptok->beta_hb);
		cout<<"beta_h "<<ptok->beta_h<<" beta_hb\t"<<ptok->beta_hb<<endl;
		tok.lookupValue("nbeta",ptok->nbeta);
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
                grid.lookupValue("nx",pgrid->nx);
                grid.lookupValue("nL",pgrid->nL);
                grid.lookupValue("nE",pgrid->nE);
                grid.lookupValue("ntheta",pgrid->ntheta);
		assert(pgrid->nx%3==1);
		assert(pgrid->nL%3==1);
		assert(pgrid->nE%3==1);
		assert(pgrid->ntheta%3==1);
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
		slowing.lookupValue("rflag",pslowing->rflag);
		slowing.lookupValue("r0",pslowing->r0);	
		slowing.lookupValue("rd",pslowing->rd);	
		Setting &set=slowing["rc"];
		int lrc =set.getLength();
		if(lrc>nrc)
		{
			cout<<"max degree is 9"<<endl;
			lrc=nrc;
		}
		for(int i=0;i<nrc;i++)
			pslowing->rc[i]=0;
		for(int i=0;i<lrc;i++)
			pslowing->rc[i] =set[i];	
			
		slowing.lookupValue("L0",pslowing->L0);	
		slowing.lookupValue("Ld",pslowing->Ld);	
		slowing.lookupValue("E0",pslowing->E0);	
		slowing.lookupValue("Ed",pslowing->Ed);	
		slowing.lookupValue("Ec",pslowing->Ec);	
		slowing.lookupValue("sigma",pslowing->sigma);
	
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
		modeset.lookupValue("n",mode->n);
		modeset.lookupValue("m",mode->m);
		modeset.lookupValue("pa",mode->pa);
		modeset.lookupValue("pb",mode->pb);
		modeset.lookupValue("delta_r",mode->delta_r);
		modeset.lookupValue("input_i",mode->input_i);
		modeset.lookupValue("input_filename",mode->mode_filename);	
		cout<<"mode filename \t"<<mode->mode_filename<<endl;
	
		modeset.lookupValue("omega_0",mode->omega_0);
		modeset.lookupValue("omega_1",mode->omega_1);
		modeset.lookupValue("omega_i",mode->omega_i);
		modeset.lookupValue("omega_n",mode->omega_n);
		modeset.lookupValue("omega_err",mode->omega_err);
		modeset.lookupValue("max_iter",mode->max_iter);
		modeset.lookupValue("max_iterg",mode->max_iterg);
		modeset.lookupValue("dw_f",mode->dw_f);
		modeset.lookupValue("zero_rhod",mode->zero_rhod);
		modeset.lookupValue("zero_iner",mode->zero_iner);
		modeset.lookupValue("xi_0",mode->xi_0);

		modeset.lookupValue("gomegar",mode->gomegar);
		modeset.lookupValue("gomegai",mode->gomegai);
		assert(mode->pa<=mode->pb);
		assert(mode->omega_n>1&&mode->omega_n<10000);
		Alloc1D(mode->omega_array,mode->omega_n);
		
		double domega=(mode->omega_1-mode->omega_0)/(mode->omega_n-1);
		std::complex<double> ti=1.0i;
		for(int iomega=0;iomega<mode->omega_n;iomega++)
		{
			mode->omega_array[iomega] = mode->omega_0 +domega*iomega +ti*mode->omega_i; 
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
		optset.lookupValue("adiabatic_i",pdwkopt->adiabatic_i);
	}	
	catch(const SettingNotFoundException &nfex)
        {
                cerr<<"No mode parameters in filename \t"<<filename<<endl;
                return(EXIT_FAILURE);
        }
	
	return 0;
}
