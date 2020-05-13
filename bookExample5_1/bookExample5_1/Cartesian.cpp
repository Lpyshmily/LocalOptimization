#include "ExternHeader.h"
#include "Constant.h"

int Derivative(double t, const double* x, double* dx, const double* dfpara)
{
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // 同伦参数
	double lam0 = dfpara[1];

	// 计算推力方向单位矢量
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	// 计算推力幅值大小
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // 对数型同伦指标下的推力

	double r = V_Norm2(rv, 3);
	for (int i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3个位置分量的导数是3个速度分量
		dx[i+3] = -muNU/(r*r*r)*rv[i] + u*TmaxNU/m*alpha[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = -u*TmaxNU/Ispg0NU;
	dx[13] = -u*TmaxNU/(m*m)*norm_lambdaV;
	return 1;
}


//利用Derivative()得到14个积分变量随时间变化的导数，通过ode45()进行给定协态初值情况下的积分，得到估计的终端状态变量
//两端时间和状态固定的燃料最优问题打靶方程
//sfpara
//[0-6],初始位置速度和质量
//[7-12],末端位置速度
//[13-14],转移时间和epsilon
//[15],输出信息标识，0表示不输出，1表示输出
//x:[0],lam0;[1-7]初始协态
//n为方程维数，fvec为残差，iflag为0时表示上方调用函数告知不必算残差
//计算成功时，返回0或大于0的值，不成功时返回小于0的值
int ComputeFvec(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if(iflag==0)
	{
		return 0;
	}
	
	double x0[14] = {0.0}, rvf[6] = {0.0};
	V_Copy(x0, sfpara, 7);
	V_Copy(&x0[7], &x[1], 7);
	double lam0 = x[0];
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	double epsi = sfpara[14];
	int outflag = sfpara[15];

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	// FILE *fid = NULL;//fopen("temp.txt","w");//如果设定有效文件路径，最后需要关闭文件
	FILE *fid = fopen("temp.txt", "w");
	int flag,NumPoint;
	double work[140]={0.0};
	flag = ode45(Derivative, x0, dfpara,  0.0, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);

	for (int i=0;i<6;++i)
		fvec[i] = x0[i] - rvf[i];
	fvec[6] = x0[13];
	fvec[7] = enorm(8, x) - 1.0;

	if(outflag>0)
	{
		fvec[0]=x0[6];
	}
	fclose(fid);
	return 0;
}

int FOP(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum)
{
	double sfpara[16] = {0.0}; // 0~6-7个初始状态值，7~13-7个终端状态值，14-同伦参数，15-ShootFunFuelOpt()的输出信息标志
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rv1, 6);
	sfpara[13] = tof;
	sfpara[14] = epsi;
	sfpara[15] = 0.0;
	
	int info, flag = 0;
	int n = 8;
	double x[8] = {0.0}, fvec[8] = {0.0}, wa[200] = {0.0};
	double xtol = 1.0e-8;
	
	int num = 0;
	while (num<MaxGuessNum)
	{
		for(int j=0;j<8;j++)
				x[j]=(double)rand()/RAND_MAX-0.5; // 随机给出打靶变量初值
		x[7]+=0.5; // 质量协态导数为负，且末端为0，因此初值必为正
		x[0]+=0.5; // lambda0，归一化乘子，人为取正
		info = hybrd1(ComputeFvec, n, x, fvec, sfpara, wa, xtol, 20, 500); // info-hybrd1()的输出标志
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			int _j = ComputeFvec(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
			if(fvec[0]>Out[0]) // 剩余质量为正，且满足打靶精度要求，停止
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 8);
				break;
			}
			sfpara[15]=0.0;
		}
		num++;
	}
	return flag;
}