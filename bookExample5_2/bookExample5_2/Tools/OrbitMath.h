#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_

#include<math.h>
#include<iostream>
#include<assert.h>
#include"..\Constant.h"


//符号函数
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

//求最大值
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}


//求最小值
template <class T>
inline T Min (T x, T y) 
{
	return (x<y)?x:y;
}

template <class T> inline T CopySign (T x, T y)
{
	return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x);
}

template <class T> void PaiXu(T* p, int N)
{

	T tmp = 0; 
	for (int i=0 ;i<N ;i++ ) 
	{ 
		for (int j=0 ;j<N-1-i ;j++ ) 
		{ 
			if (p[j] > p[j+1]) 
			{ 
				tmp = p[j]; 
				p[j] = p[j+1]; 
				p[j+1] = tmp;
			}
		} 
	}
}

//缔合勒让德多项式
double LegendreP(int BigL, int SmallM, double Variable);

//归一化的勒让德多项式值
double NormLegendreP(int BigL, int SmallM, double Variable);

//归一化的勒让德多项式值Pnm(x)组成的二维数组,行i=0,1,2,....,n,列j=0,1,2,...,m.当j>i时,Pnm(x)=0
void NormLPArray(int n, int m, double** NormPNM, double x);



//取最近接x的整数,返回double型.
inline double ANINT(double x)
{
	double left=fmod(x, 1.0);
	if(fabs(left)<0.5)
		return x-left;
	else if(left>=0.5)
		return x-left+1;
	else
		return x-left-1;
}

//取最近接x的整数,返回int型.
inline int NINT(double x)
{
	double left=ANINT(x);
	return (int)left;
}

inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

inline double acosh(double x)
{	
	assert(x>=1.0);
	return log(x+sqrt(x*x-1.0));
}

inline double randn(double mu,double sigma)
{
	const int n=6;
	assert(sigma>=0);
	if(sigma==0) return mu;
	double sum=0;
	for(int i=0;i<2*n;i++)
		sum+=(double)rand()/(double)RAND_MAX;//g_rg.Random();
	return (sum-n)*sigma+mu;
}
//Lagrange插值函数
double LagInterp(double* x, double* y, int order, double t);


#endif