#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_

#include <assert.h>
#include <cmath>

// ���ź���
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

// ���������нϴ��
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y) ? x : y;
}


// ���������н�С��
template <class T>
inline T Min (T x, T y) 
{
	return (x<y) ? x : y;
}

// ȡ�����x������,����double��.
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

// ȡ�����x������,����int��.
inline int NINT(double x)
{
	double left=ANINT(x);
	return (int)left;
}

// ��˫������ֵ
inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

// ��˫������ֵ
inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

#endif