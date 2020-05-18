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
// 根据经典轨道根数求地心惯性直角坐标系下的位置和速度分量
void coe2rv(int& flag, double* rv, const double* coe, double mu=3.98600441800e+14);
// E2f 根据偏近点角和偏心率求真近点角
double E2f(int& flag, double E, double e);
// E2M 根据偏近点角和偏心率求平近点角
double E2M(int& flag, double E, double e);
// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu=3.98600441800e+14, int MaxIter=100, double epsilon=1.0e-14);
// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu=3.98600441800e+14);
// f2E 根据真近点角和偏心率求偏近点角
double f2E(int& flag, double f, double e);
// Inertia2LVLH 根据主星的绝对位置R0和速度V0及从星的绝对位置R1和速度V1求从星在主星轨道坐标系中的相对位置和速度
void Inertia2LVLH(int& flag, double* rv, const double* RV0, const double* RV1);
// LVLH2Inertia 根据主星的绝对位置R0和速度V0及从星相对主星轨道坐标系中的位置r和速度v求从星的绝对位置和速度
void LVLH2Inertia(int& flag, double* RV1, const double* RV0, const double* rv);
// M2E 根据平近点角和偏心率求偏近点角
double M2E(int& flag, double M, double e, int MaxIter=100, double epsilon=1.0e-14);
// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
void rv2coe(int& flag, double* coe, const double* RV, double mu=3.98600441800e+14);
//Lambert问题求解
void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2, 
			 double tf, const double* unith, int& flag, double mu=3.98600441800e+14, int way=0, int N=0, int branch=0, 
			 int Maxiter=60, double tol=1.0e-12);
//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int&flag, double* rv1, const double* rv0, double t0, double t1, double mu=3.98600441800e+14);
//根据初始状态rv0和末端状态rv1及转移时间t，求特征速度最小的双脉冲转移轨道。
//返回：dv0,初始脉冲矢量;dv1，末端脉冲矢量;Mdv0，dv0的幅值；Mdv1,dv1的幅值；N，最优圈次;branch，最优分枝;flag，1表示计算成功
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
			  const double* rv0, const double* rv1, double t, double mu=3.98600441800e+14);

// dee_dt 求天球轨道根数ee在沿径向、横向、法向加速度Ad作用下随时间的变化率
void dee_dt(int&flag, double* dee, const double* Ad, const double* ee, double mu=3.98600441800e+14);
// dmee_dt 求修正天球轨道根数mee在沿径向、横向、法向加速度Ad作用下随时间的变化率
void dmee_dt(int& flag, double* dmee, const double* Ad, const double* mee, int orbtype=1, double mu=3.98600441800e+14);
// ee2coe 根据天球根数求经典轨道根数
void ee2coe(int&flag, double* coe, const double* ee, double mu=3.98600441800e+14);
// ee2rv 根据天球轨道根数求地心惯性直角坐标系下的位置和速度分量
void ee2rv(int&flag, double* rv, const double* ee, double mu=3.98600441800e+14);
// coe2ee 根据经典轨道根数求天球轨道根数，对轨道倾角180度奇异
void coe2ee(int&flag, double* ee, const double* coe, double mu=3.98600441800e+14);
// coe2mee 根据经典轨道根数求修正天球轨道根数
void coe2mee(int& flag, double* mee, int& orbtype, const double* coe, double mu=3.98600441800e+14);
// mee2coe 根据修正天球轨道根数求经典轨道根数
void mee2coe(int&flag, double* coe, const double* mee, int orbtype=1, double mu=3.98600441800e+14);
// mee2rv 根据修正天球轨道根数求地心惯性直角坐标系下的位置和速度分量
// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求天球轨道根数，对轨道倾角180度时奇异
void rv2ee(int&flag, double* ee, const double* RV, double mu=3.98600441800e+14);
// rv2mee 根据地心惯性直角坐标系下的位置和速度分量求修正天球轨道根数
void rv2mee(int&flag, double* mee, int& orbtype, const double* RV, double mu=3.98600441800e+14);

#endif