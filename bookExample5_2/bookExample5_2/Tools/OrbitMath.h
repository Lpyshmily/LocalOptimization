#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_

#include<math.h>
#include<iostream>
#include<assert.h>
#include"..\Constant.h"


//���ź���
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

//�����ֵ
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}


//����Сֵ
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

//�޺����õ¶���ʽ
double LegendreP(int BigL, int SmallM, double Variable);

//��һ�������õ¶���ʽֵ
double NormLegendreP(int BigL, int SmallM, double Variable);

//��һ�������õ¶���ʽֵPnm(x)��ɵĶ�ά����,��i=0,1,2,....,n,��j=0,1,2,...,m.��j>iʱ,Pnm(x)=0
void NormLPArray(int n, int m, double** NormPNM, double x);



//ȡ�����x������,����double��.
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

//ȡ�����x������,����int��.
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
//Lagrange��ֵ����
double LagInterp(double* x, double* y, int order, double t);


#endif