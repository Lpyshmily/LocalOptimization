#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_

#include <math.h>
#include "..\Constant.h"
#include "OrbitMath.h"
#include "VecMat.h"

// E2f ����ƫ����Ǻ�ƫ������������
double E2f(int& flag, double E, double e);
// E2M ����ƫ����Ǻ�ƫ������ƽ�����
double E2M(int& flag, double E, double e);

// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu=3.98600441800e+14, int MaxIter=100, double epsilon=1.0e-14);
// f0ft2dt ���ݳ�ʼ�����Ǻ��������������ݻ�ʱ��
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu=3.98600441800e+14);

// f2E ���������Ǻ�ƫ������ƫ�����
double f2E(int& flag, double f, double e);

// M2E ����ƽ����Ǻ�ƫ������ƫ����ǣ��⿪���շ���
double M2E(int& flag, double M, double e, int MaxIter=100, double epsilon=1.0e-14);

// rv2ee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ������������������Թ�����180��ʱ����
void rv2ee(int&flag, double* ee, const double* RV, double mu=3.98600441800e+14);
#endif