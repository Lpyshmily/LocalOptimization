#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_

#include<math.h>
#include<assert.h>
#include<iostream>
#include"..\Constant.h"
#include<iomanip>
#include"OrbitMath.h"
#include"VecMat.h"
//namespace OrbitFun{
// ���ݾ�������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
void coe2rv(int& flag, double* rv, const double* coe, double mu=3.98600441800e+14);
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
// Inertia2LVLH �������ǵľ���λ��R0���ٶ�V0�����ǵľ���λ��R1���ٶ�V1����������ǹ������ϵ�е����λ�ú��ٶ�
void Inertia2LVLH(int& flag, double* rv, const double* RV0, const double* RV1);
// LVLH2Inertia �������ǵľ���λ��R0���ٶ�V0������������ǹ������ϵ�е�λ��r���ٶ�v����ǵľ���λ�ú��ٶ�
void LVLH2Inertia(int& flag, double* RV1, const double* RV0, const double* rv);
// M2E ����ƽ����Ǻ�ƫ������ƫ�����
double M2E(int& flag, double M, double e, int MaxIter=100, double epsilon=1.0e-14);
// rv2coe ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ����󾭵�������
void rv2coe(int& flag, double* coe, const double* RV, double mu=3.98600441800e+14);
//Lambert�������
void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2, 
			 double tf, const double* unith, int& flag, double mu=3.98600441800e+14, int way=0, int N=0, int branch=0, 
			 int Maxiter=60, double tol=1.0e-12);
//���ݳ�ʼʱ��t0��״̬rv0��ĩ��ʱ��t1��״̬rv1���������ƽ���������ɹ�,flag����1
void rv02rvf(int&flag, double* rv1, const double* rv0, double t0, double t1, double mu=3.98600441800e+14);
//���ݳ�ʼ״̬rv0��ĩ��״̬rv1��ת��ʱ��t���������ٶ���С��˫����ת�ƹ����
//���أ�dv0,��ʼ����ʸ��;dv1��ĩ������ʸ��;Mdv0��dv0�ķ�ֵ��Mdv1,dv1�ķ�ֵ��N������Ȧ��;branch�����ŷ�֦;flag��1��ʾ����ɹ�
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
			  const double* rv0, const double* rv1, double t, double mu=3.98600441800e+14);

// dee_dt ������������ee���ؾ��򡢺��򡢷�����ٶ�Ad��������ʱ��ı仯��
void dee_dt(int&flag, double* dee, const double* Ad, const double* ee, double mu=3.98600441800e+14);
// dmee_dt ����������������mee���ؾ��򡢺��򡢷�����ٶ�Ad��������ʱ��ı仯��
void dmee_dt(int& flag, double* dmee, const double* Ad, const double* mee, int orbtype=1, double mu=3.98600441800e+14);
// ee2coe ������������󾭵�������
void ee2coe(int&flag, double* coe, const double* ee, double mu=3.98600441800e+14);
// ee2rv ������������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
void ee2rv(int&flag, double* rv, const double* ee, double mu=3.98600441800e+14);
// coe2ee ���ݾ������������������������Թ�����180������
void coe2ee(int&flag, double* ee, const double* coe, double mu=3.98600441800e+14);
// coe2mee ���ݾ�������������������������
void coe2mee(int& flag, double* mee, int& orbtype, const double* coe, double mu=3.98600441800e+14);
// mee2coe �������������������󾭵�������
void mee2coe(int&flag, double* coe, const double* mee, int orbtype=1, double mu=3.98600441800e+14);
// mee2rv ����������������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
// rv2ee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ������������������Թ�����180��ʱ����
void rv2ee(int&flag, double* ee, const double* RV, double mu=3.98600441800e+14);
// rv2mee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ�������������������
void rv2mee(int&flag, double* mee, int& orbtype, const double* RV, double mu=3.98600441800e+14);

#endif