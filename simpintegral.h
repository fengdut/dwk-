#ifndef SIMPINTEGRAL_H
#define SIMPINTEGRAL_H

//Simpson's 3/8 rule
//ny mush equiv to 3*n+1
//double simpintegral(double *y,int ny,double dx);

template <class datatype>
datatype simpintegral(datatype *y, int ny, double dx)
{
        datatype sum=0;
        for(int i=0;i<ny;i++)
        {
                sum +=3.0*y[i];
        }
        for(int i=3;i<ny-1;i=i+3)
        {
                sum -= y[i];
        }
        sum =sum -2.0*(y[0]+y[ny-1]);
        sum =dx*3.0*sum*0.125;
        return sum;

}

#endif

