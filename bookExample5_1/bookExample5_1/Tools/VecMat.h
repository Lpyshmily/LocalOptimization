#ifndef _VECMAT_H_
#define _VECMAT_H_
#include <iostream>
#include <math.h>
#include <assert.h>
using namespace std;
// �ͷ��ڴ�
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B = NULL;
}


// ������A��ֵ����B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) B[I_]=A[I_];
}

// ������ֵ���θ���B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0] = x;
	B[1] = y;
	B[2] = z;
}

// ������ֵ���θ���B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0] = x;
	B[1] = y;
	B[2] = z;
	B[3] = vx;
	B[4] = vy;
	B[5] = vz;
}


// �ж�����B������Aÿ��Ԫ���Ƿ���ȣ����ʱ������
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) if(B[I_] != A[I_]) return false;
	return true;
}


// ȡ����������B��ֵ��ʹB[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;++I_) B[I_]=-A[I_];
}


// �ӷ���������C��ֵ��ʹC[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] + B[I_];
}

// �ӷ���������C��ֵ��ʹC[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A + B[I_];	
}


// ������������C��ֵ��ʹC[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] - B[I_];
}

// ������������C��ֵ��ʹC[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] - B;	
}


// �˷���������C��ֵ��ʹC[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] * B[I_];
}

// �˷���������C��ֵ��ʹC[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] * B;	
}


// ������������C��ֵ��ʹC[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] / B[I_];
}

// ������������C��ֵ��ʹC[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;++I_) C[I_] = A[I_] / B;	
}


// ���ڻ�������
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;++I_) result += A[I_] * B[I_];
	return result;
}

// �����C=AXB,������V_Cross(B,B,A)��V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0] = A[1] * B[2] - A[2] * B[1];
	C[1] = A[2] * B[0] - A[0] * B[2];
	C[2] = A[0] * B[1] - A[1] * B[0];
}

// ȡ����ֵ��A[i]=|B[i]|
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

// ��1-����������
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

// ��2-����������
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result = V_Dot(B,B,N);
	return sqrt(result);
}

// ������-����������
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

// �����Ԫ�أ�������
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;++I_) if(B[I_]>result) result = B[I_];
	return result;
}

// �����Ԫ�أ������أ��Ҹı�������
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;++I_) if(B[I_]>maximal) {index = I_; maximal = B[I_];}
	return maximal;
}

// ����СԪ�أ�������
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;++I_) if(B[I_]<result) result = B[I_];
	return result;
}

// ����СԪ�أ������أ��Ҹı�������
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;++I_) if(B[I_]<minimal) {index = I_; minimal = B[I_];}
	return minimal;
}

// ����������
template<class T> inline void V_Input(istream  &input, T* Vec, int N)
{
    for(int I_=0;I_<N;++I_)	input >> Vec[I_];
}

// ���������
template<class T> inline void V_Output(ostream &output, const T* Vec, int N)
{
    for(int I_=0;I_<N;++I_) output << Vec[I_] << endl;   
}
#endif