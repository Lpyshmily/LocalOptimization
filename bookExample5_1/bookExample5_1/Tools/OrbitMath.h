#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_

#include <assert.h>
#include <cmath>

// 符号函数
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

// 求两个数中较大的
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y) ? x : y;
}


// 求两个数中较小的
template <class T>
inline T Min (T x, T y) 
{
	return (x<y) ? x : y;
}

// 取最近接x的整数,返回double型.
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

// 取最近接x的整数,返回int型.
inline int NINT(double x)
{
	double left=ANINT(x);
	return (int)left;
}

// 求双曲正切值
inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

// 求双曲正弦值
inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

#endif