#ifndef _VECMAT_H_
#define _VECMAT_H_
#include <iostream>
#include <math.h>
#include <assert.h>
using namespace std;
// 释放内存
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B = NULL;
}


// 将向量A的值赋给B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) B[I_]=A[I_];
}

// 将三个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0] = x;
	B[1] = y;
	B[2] = z;
}

// 将六个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0] = x;
	B[1] = y;
	B[2] = z;
	B[3] = vx;
	B[4] = vy;
	B[5] = vz;
}


// 判断向量B与向量A每个元素是否相等，相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) if(B[I_] != A[I_]) return false;
	return true;
}


// 取反，给向量B赋值，使B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) B[I_]=-A[I_];
}


// 加法，给向量C赋值，使C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] + B[I_];
}

// 加法，给向量C赋值，使C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A + B[I_];	
}


// 减法，给向量C赋值，使C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] - B[I_];
}

// 减法，给向量C赋值，使C[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] - B;	
}


// 乘法，给向量C赋值，使C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] * B[I_];
}

// 乘法，给向量C赋值，使C[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] * B;	
}


// 除法，给向量C赋值，使C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] / B[I_];
}

// 除法，给向量C赋值，使C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] / B;	
}


// 求内积并返回
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;++I_) result += A[I_] * B[I_];
	return result;
}

// 求外积C=AXB,不能用V_Cross(B,B,A)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0] = A[1] * B[2] - A[2] * B[1];
	C[1] = A[2] * B[0] - A[0] * B[2];
	C[2] = A[0] * B[1] - A[1] * B[0];
}

// 取绝对值，A[i]=|B[i]|
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_)
	{
		if(B[I_] >= 0)
			A[I_] = B[I_];
		else
			A[I_] = -B[I_];
	}
}

// 求1-范数并返回
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result = 0;
	for(int I_=0;I_<N;++I_)
	{
		if(B[I_] >= 0)
			result += B[I_];
		else
			result -= B[I_];
	}
	return result;
}

// 求2-范数并返回
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result = V_Dot(B,B,N);
	return sqrt(result);
}

// 求无穷-范数并返回
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result = 0;
	for(int I_=0;I_<N;++I_) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}

// 求最大元素，并返回
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;++I_) if(B[I_]>result) result = B[I_];
	return result;
}

// 求最大元素，并返回，且改变索引数
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;++I_) if(B[I_]>maximal) {index = I_; maximal = B[I_];}
	return maximal;
}

// 求最小元素，并返回
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;++I_) if(B[I_]<result) result = B[I_];
	return result;
}

// 求最小元素，并返回，且改变索引数
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;++I_) if(B[I_]<minimal) {index = I_; minimal = B[I_];}
	return minimal;
}

// 输入流函数
template<class T> inline void V_Input(istream  &input, T* Vec, int N)
{
    for(int I_=0;I_<N;++I_)	input >> Vec[I_];
}

// 输出流函数
template<class T> inline void V_Output(ostream &output, const T* Vec, int N)
{
    for(int I_=0;I_<N;++I_) output << Vec[I_] << endl;   
}
#endif