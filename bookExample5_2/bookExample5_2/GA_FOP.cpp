#include "ExternHeader.h"
#include "GA_FOP.h"

// 根据14个积分变量x，计算其对应的导数dx
// dfpara:[0],epsi;[1],lam0
// 第一个参数t只是为了ode45调用，在程序中不起作用
int GA_derivative(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // 同伦参数
	double lam0 = dfpara[1];

	// 计算推力方向单位矢量和推力幅值大小
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // 对数型同伦指标下的推力

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
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

double GA_Hamilton(const double* x, double lam0, double epsi)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// 计算推力方向单位矢量和推力幅值大小
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // 对数型同伦指标下的推力
	// 计算状态变量的导数
	double dfpara[2];
	dfpara[0] = epsi;
	dfpara[1] = lam0;
	double dx[14] = {0.0};
	int flag = GA_derivative(0, x, dx, dfpara);

	double H = 0.0;
	H = V_Dot(costate, dx, 7);
	H = H + (lam0*TmaxNU/Ispg0NU*( u-epsi*log(u*(1.0-u)) ));
	return H;
}

//利用GA_derivative()得到14个积分变量随时间变化的导数，通过ode45()进行给定协态初值情况下的积分，得到估计的终端状态变量
//两端时间和状态固定的燃料最优问题打靶方程
//sfpara
//[0-6],初始位置速度和质量
//[7-12],末端位置速度
//[13-14],转移时间和epsilon
//[15],输出信息标识，0表示不输出，1表示输出
//[16-21],火星在MJD_MARS时刻的位置速度
//x:[0],lam0;[1-7]初始协态;[8-11]等式约束乘子;[12]不等式约束乘子;[13-15]引力辅助速度增量dvg;[16]引力辅助时刻
//n为方程维数，fvec为残差，iflag为0时表示上方调用函数告知不必算残差
//计算成功时，返回0或大于0的值，不成功时返回小于0的值
int GA_fvec(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if(iflag==0)
	{
		return 0;
	}
	
	double x0[14], rvf[6];
	V_Copy(x0, sfpara, 7);
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	double epsi = sfpara[14];
	int outflag = NINT(sfpara[15]);
	// 获取tm时刻火星位置
	double rv_mars[6];
	V_Copy(rv_mars, &sfpara[16], 6);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7);
	double chi[4], kappa, dvg[3], tm;
	V_Copy(chi, &x[8], 4);
	kappa = x[12];
	V_Copy(dvg, &x[13], 3);
	tm = x[16];
	



	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	// FILE *fid = NULL;//fopen("temp.txt","w");//如果设定有效文件路径，最后需要关闭文件
	FILE *fid1 = fopen("temp1.txt", "w");
	FILE *fid2 = fopen("temp2.txt", "w");
	int flag,NumPoint;
	// double work[140]={0.0};
	double work[200]={0.0};
	flag = ode45(GA_derivative, x0, dfpara,  0.0, tm*TOF*86400/TUnit, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid1);
	// 计算一部分偏差
	// 引力辅助位置约束
	double rv_tm[6];
	rv02rvf(flag, rv_tm, rv_mars, MJD_MARS*86400/TUnit, (MJD0+tm*TOF)*86400/TUnit, muNU);
	for (int i=0;i<3;++i)
		fvec[7+i] = x0[i] - rv_tm[i];
	// normvinfin==normvinfout
	double vinfin[3], vinfout[3], normvinfin, normvinfout, unit_in[3], unit_out[3];
	V_Minus(vinfin, &x0[3], &rv_tm[3], 3);
	V_Add(vinfout, vinfin, dvg, 3);
	normvinfin = V_Norm2(vinfin, 3);
	normvinfout = V_Norm2(vinfout, 3);
	V_Divid(unit_in, vinfin, normvinfin, 3);
	V_Divid(unit_out, vinfout, normvinfout, 3);
	fvec[10] = normvinfin - normvinfout;
	// 互补松弛条件
	double theta = acos(V_Dot(unit_in, unit_out, 3));
	double rp = mpp/(normvinfin*normvinfout) * (1/sin(theta/2) - 1);
	// fvec[11] = kappa*(1 - rp/rmin);
	fvec[11] = 1 - rp/rmin;

	double A[3], B[3], C[3], temp, c;
	temp = 1/( 4*sin(theta/2)*sin(theta/2) * (1-sin(theta/2)) );
	for (int i=0;i<3;++i)
	{
		A[i] = rp/normvinfin*(temp*(unit_out[i] - cos(theta)*unit_in[i]) - unit_in[i]);
		B[i] = rp/normvinfout*(temp*(unit_in[i] - cos(theta)*unit_out[i]) - unit_out[i]);
		C[i] = rp*( -temp*(1/normvinfout*(unit_in[i] - cos(theta)*unit_out[i]) + 1/normvinfin*(unit_out[i] - cos(theta)*unit_in[i])) + unit_out[i]/normvinfout + unit_in[i]/normvinfin );
	}
	double am[3];
	double r_tm = V_Norm2(rv_tm, 3);
	for (int i=0;i<3;++i)
		am[i] = -muNU/(r_tm*r_tm*r_tm)*rv_tm[i];
	c = V_Dot(C, am, 3);

	// tm-时刻速度协态
	for (int i=0;i<3;++i)
		fvec[12+i] = x0[10+i] - chi[3]*unit_in[i] + 1/rmin*kappa*A[i];

	double H1 =  GA_Hamilton(x0, lam0, epsi);
	// 计算新的状态变量和协态变量
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // 速度
		x0[7+i] = x0[7+i] - chi[i]; // 位置协态
		x0[10+i] = chi[3]*unit_out[i] + 1/rmin*kappa*B[i]; // 速度协态
	}
	double H2 = GA_Hamilton(x0, lam0, epsi);
	// 静态条件偏差
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - 1/rmin*kappa*c;


	// 第二阶段的积分
	// flag = ode45(GA_derivative, x0, dfpara,  tm*TOF*86400/TUnit,  tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	flag = ode45(GA_derivative, x0, dfpara, 0.0, tf-tm*TOF*86400/TUnit, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);

	for (int i=0;i<6;++i)
		fvec[i] = x0[i] - rvf[i];
	fvec[6] = x0[13];
	fvec[16] = enorm(13, x) - 1.0;

	if(outflag>0)
	{
		fvec[0]=x0[6];
	}
	fclose(fid1);
	fclose(fid2);
	return 0;
}

int GA_FOP(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle)
{
	double sfpara[22] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rv1, 6);
	sfpara[13] = tof;
	sfpara[14] = epsi;
	sfpara[15] = 0.0;
	V_Copy(&sfpara[16], rv_middle, 6);

	
	int info, flag = 0;
	const int n = 17;
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[500] = {0.0};
	double xtol = 1.0e-8;
	
	int num = 0;
	double amp, phi, delta;
	while (num<MaxGuessNum)
	{
		for(int j=0;j<17;j++)
				x[j]=(double)rand()/RAND_MAX-0.5; // 随机给出打靶变量初值
		x[0] += 0.5; // lambda0，归一化乘子，人为取正
		x[7] += 0.5; // 质量协态导数为负，且末端为0，因此初值必为正
		amp = x[13];
		phi = (x[14] - 0.5)*DPI;
		delta = x[15]*D2PI;
		x[13] = sqrt(mpp/rmin)*amp*cos(phi)*cos(delta);
		x[14] = sqrt(mpp/rmin)*amp*cos(phi)*sin(delta);
		x[15] = sqrt(mpp/rmin)*amp*sin(phi);
		x[16] += 0.5; // 引力辅助时间也要是正的
		info = hybrd1(GA_fvec, n, x, fvec, sfpara, wa, xtol, 20, 500); // info-hybrd1()的输出标志
		std::cout << n << std::endl;
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			int _j = GA_fvec(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
			if(fvec[0]>Out[0]) // 剩余质量为正，且满足打靶精度要求，停止
			{
				flag=1;
				Out[0]=fvec[0];
				V_Copy(&Out[1], x, 17);
				break;
			}
			sfpara[15]=0.0;
		}
		num++;
	}
	return flag;
}