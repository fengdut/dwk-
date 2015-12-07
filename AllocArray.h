#ifndef _ALLOCARRAY_H
#define _ALLOCARRAY_H

#include <memory.h>
#include<assert.h>
template<class datatype>
void Alloc2D(datatype **& pdata,int n,int m)
{
	pdata	=	new datatype*[n];
	pdata[0]=	new datatype[n*m];
	memset(pdata[0],0,sizeof(datatype)*n*m);

	for(int i=1;i<n;i++)
	{
		pdata[i]	=	pdata[i-1]+m;
	}
};

template<class datatype>
void Free2D(datatype **& pdata)
{
	delete[] pdata[0];
	delete[] pdata;
};


template<class datatype>
void Alloc1D(datatype *& pdata,size_t n)
{
	pdata= new datatype[n];
}

template<class datatype>
void Free1D(datatype *& pdata)
{
	delete[] pdata;
	pdata=0;
}


template<class datatype>
void Alloc3D(datatype ***& pdata,size_t m, size_t n, size_t z)
{
	assert(m*n*z<100000000);
	pdata	=	new datatype**[m];

	pdata[0]=	new datatype*[n*m];
	for(size_t i=1;i<m;i++)
		pdata[i]=pdata[i-1]+n;

	pdata[0][0]=new datatype[m*n*z];
	for(size_t i=1;i<m;i++)
		pdata[i][0]=pdata[i-1][0]+n*z;

	for(size_t i=0;i<m;i++)
		for(size_t j=1;j<n;j++)
		{
			pdata[i][j]=pdata[i][j-1]+z;
		}
}

template <class datatype>
void Free3D(datatype ***&pdata)
{
	delete[] pdata[0][0];
	delete[] pdata[0];
	delete[] pdata;
}

#endif
