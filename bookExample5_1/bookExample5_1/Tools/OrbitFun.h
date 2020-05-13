#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_

#include <math.h>
#include "..\Constant.h"
#include "OrbitMath.h"
#include "VecMat.h"

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

// M2E 根据平近点角和偏心率求偏近点角，解开普勒方程
double M2E(int& flag, double M, double e, int MaxIter=100, double epsilon=1.0e-14);

// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求天球轨道根数，对轨道倾角180度时奇异
void rv2ee(int&flag, double* ee, const double* RV, double mu=3.98600441800e+14);
#endif