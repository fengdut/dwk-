#include"stdio.h"
#include"dwk.h"
#include<complex>
#include"math.h"
void fundwk(double omegai,double *x, double *f, double *A,complex<double> &dwk);
void fundwk2(double *x, double *f, double *A, complex<double> &dw);

int newton_iter(double omegai,double& omegar,double& betah,complex<double> &dwk)
{
	int maxit=100;
	double eps=1e-12;
	double x[2],f[2],A[4];
	x[0]	=	betah;
	x[1]	=	omegar;
	double domega=0.005;
	double dbeta=0.001;
	double delta[2]={0};
	double B[2]={0};
	cout<<scientific;
	int i=0;
	int j=1;
	cout<<"init omegar:\t"<<omegar<<"\t betah:\t"<<betah<<endl;
	for(i=0;i<maxit;i++)
	{
		fundwk(omegai,x,f,A,dwk);		
		if(abs(f[0])<eps&&abs(f[1])<eps)
		{
			if(x[1]>0&&x[1]<0.1&&x[0]>0&&x[0]<1)
			{
			cout<<"f[0] \t"<<f[0]<<"\t f[1] \t"<<f[1]<<endl;
			cout<<"newton iter converged"<<endl;
			omegar=x[1];
			betah=x[0];
			break;
			}
			else
			{
				x[1]=omegar+domega*j;
				x[0]=betah+dbeta*j;
				j=j+1;
				cout<<"******* reset init point*******"<<endl;
				cout<<"omega0:\t"<<x[1]<<"\t betah0:"<<x[0]<<endl;
			}
		}
		else
		{
			B[0]	=-f[0];
			B[1]	=-f[1];
			delta[1]	=	(B[1]*A[0]-A[2]*B[0])/(A[0]*A[3]-A[1]*A[2]);
			delta[0]	=	(B[0]-A[1]*delta[1])/A[0];
			x[0]	+=delta[0];
			x[1]	+=delta[1];
			cout<<"***iter No. "<<i<<"\tbeta\t"<<x[0]<<"\tomega \t"<<x[1]<<endl;
			cout<<"f: \t"<<f[0]<<"\t"<<f[1]<<endl;
	//		cout<<"delta: \t"<<delta[0]<<"\t"<<delta[1]<<endl;
		}

	}	
	if(i==maxit)
	{
		cout<<"!!!!!!!!!newton iter not converged!!!!!!"<<endl;
		return -1;
	}
	cout<<fixed;
	return 0;
}

int newton_iter2(double &omegai,double &omegar,complex<double> &dw)
{
	int maxit	=50;
	double eps=1e-15;
	double x[2],f[2],A[4];
	x[0]=omegar;
	x[1]=omegai;
	double delta[2]={0};
	double B[2]={0};
	cout<<scientific;
	int i=0;
	cout<<"init **** omegar\t"<<x[0]<<"\tomegai\t"<<x[1]<<endl;
	for(i=0;i<maxit;i++)
	{
		fundwk2(x,f,A,dw);
		if(abs(f[0])<eps&&abs(f[1])<eps)
		{
			omegai=x[1];
			omegar=x[0];
			cout<<"omegar\t"<<x[0]<<endl;
			cout<<"omegai\t"<<x[1]<<endl;
			return 0;
		}
		B[0]    =-f[0];
		B[1]    =-f[1];
		 delta[1]        =       (B[1]*A[0]-A[2]*B[0])/(A[0]*A[3]-A[1]*A[2]);
		 delta[0]        =       (B[0]-A[1]*delta[1])/A[0];
		 x[0]    +=delta[0];
		 x[1]    +=delta[1];
		 cout<<"iter No."<<i<<endl;
		 cout<<"omegar\t"<<x[0]<<"\tomegai\t"<<x[1]<<endl;
		 cout<<"f: \t"<<f[0]<<"\t"<<f[1]<<endl;
	}
	return -1;

}

void fundwk(double omegai,double *x, double *f, double *A,complex<double> &dwk)
{
	complex<double> ti;
	ti.real(0.0);
	ti.imag(1.0);
	double small=0.000001;
	
	complex<double> omega=x[1] +omegai*ti;
//	cout<<"omega"<<omega<<endl;
	dwk=tdwk_omega(omega);
	complex<double> dwk0=tdwk_omega(omega-small);
	complex<double> dwk1=tdwk_omega(omega+small);
	

	f[0]	=omegai	+	x[0]*dwk.real();
	f[1]	=x[1]	-	x[0]*dwk.imag();

	A[0]	=dwk.real();
	A[1]	=x[0]*(dwk1.real()-dwk0.real())/(2*small);

	A[2]	=-1.0*dwk.imag();
	A[3]	=1-x[0]*(dwk1.imag()-dwk0.imag())/(2*small);
//	cout<<"x \t"<<x[0]<<"\t"<<x[1]<<endl;	
//	cout<<"A\t"<<A[0]<<"\t"<<A[1]<<"\t"<<A[2]<<"\t"<<A[3]<<endl;
//	cout<<"f\t"<<f[0]<<"\t"<<f[1]<<endl;
	
}

void fundwk2(double *x, double *f, double *A,complex<double>& dw)
{
	complex<double> ti;
        ti.real(0.0);
        ti.imag(1.0);
        double small=0.000001;
	complex<double> omega=x[0] +0.0003*ti;

	dw=tdw_omega(omega);

	f[0]	=	x[1]+dw.real();
	f[1]	=	x[0]-dw.imag();

	complex<double>	dw0	=tdw_omega(omega-small);
	complex<double> dw1	=tdw_omega(omega+small);
	A[0]	=(dw1.real()-dw0.real())/(2*small);

//	complex<double> dwa     =tdw_omega(omega-small*ti);
//	complex<double> dwb     =tdw_omega(omega+small*ti);
//	A[1]	=1	+	(dwb.real()-dwa.real())/(2*small);
	A[1]	=1;

	A[2]	=1	- (dw1.imag()-dw0.imag())/(2*small);
//	A[3]	=	-(dwb.imag()-dwa.imag())/(2*small);
	A[3]	=	0;
}


