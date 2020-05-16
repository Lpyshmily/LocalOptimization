#ifndef _MODELFUN_H_
#define _MODELFUN_H_

// 已知春分点轨道根数ee及第6个角度量初值L0和末值Lt，求从L0变化到Lt所需的时间，mus为日心引力系数，按开普勒轨道运行
double L0Lt2dt(double L0, double Lt, const double* ee, double mus);

// 已知春分点轨道根数ee及第6个角度量初值L0，求经过时间dt后的Lt，mus为日心引力系数，按开普勒轨道运行
double L0dt2Lt(double L0, double dt, const double* ee, double mus);

//求春分点轨道根数ee描述下动力学方程中的M矩阵(6*3)，D矩阵第6个分量值D6，
//M矩阵及D6对ee的偏导数三阶矩阵dM(6*3*6)和一阶数组dD6(6)
//M矩阵及D6对ee的二阶偏导数三阶矩阵ddM(6*3*6*6)和二阶数组ddD6(6*6)
//step:1,只算M, D6; 2, 算到dM, dD6;3, 算到ddM, ddD6
void MatMD(double * M, double& D6, double * dM, double* dD6, double* ddM, double* ddD6, const double* ee, int step);

//燃料最优问题动力学微分方程，包括状态量和协态量，共14维
//t为时间,x为当前状态和协态,dx为随时间的变化率, dfpara是所需的参数集
//计算成功时，返回值为1，否则返回小于1的整数值
int DynEquFuelOpt(double t, const double* x, double* dx, const double* dfpara);

//两端时间和状态固定的燃料最优问题打靶方程
//sfpara
//[0-6],初始春分点轨根和质量
//[7-12],末端春分点轨根
//[13-14],转移时间和同伦参数epsilon
//[15],输出信息标识
//x:[0],lam0;[1-7]初始协态
//n为方程维数，fvec为残差，iflag为0时表示上方调用函数告知不必算残差
//计算成功时，返回0或大于0的值，不成功时返回小于0的值
int ShootFunFuelOpt(int n, const double *x, double *fvec, int iflag, const double* sfpara);

//反复随机猜初值，通过ShootFunFuelOpt()计算给定初值条件下的积分结果，再通过hybrid()进行打靶方程求解
//求两端固定的燃料最优问题
//Out为输出结果，分别为剩余质量，8个乘子初值（第一个为lamd0)，共9维
//ee0为出发的春分点轨道根数，ee1为到达的春分点轨道根数，m0为出发质量，ft为飞行时间，
//epsi为同伦参数的值（直接取一个较小的值，如0.001，求得接近邦邦控制的结果），MaxGuessNum为最大猜测次数
//计算成功时，返回1，否则返回0
int Fixed2PFuelOpt(double* Out, const double* ee0, const double* ee1, double m0, double ft, double epsi, int MaxGuessNum);

void OutputAmplitude(const double* Out, const double* ee0, double m0, double epsi, const double* TimeNode, int TimeNodeNum);

#endif