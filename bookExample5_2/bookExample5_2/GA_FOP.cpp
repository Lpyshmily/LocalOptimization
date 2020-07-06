#include "ExternHeader.h"
#include "GA_FOP.h"
#include "integrate.h"

// 对数同伦
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

// bang-bang控制
// 根据14个积分变量x，计算其对应的导数dx
int GA_derivative_type1(double t, const double* x, double* dx, const double* dfpara)
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
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	if (rou > 0)
		u = 0.0;
	else
		u = 1.0;
	

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

// 二次同伦
int GA_derivative_type2(double t, const double* x, double* dx, const double* dfpara)
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
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/2.0/epsi;

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

// 对数同伦
// 计算哈密顿函数的值
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

// bang-bang控制
// 计算哈密顿函数的值
double GA_Hamilton_type1(const double* x, double lam0, double epsi)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// 计算推力方向单位矢量和推力幅值大小
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	if (rou > 0.0)
		u = 0.0;
	else
		u = 1.0;
	// 计算状态变量的导数
	double dfpara[2];
	dfpara[0] = epsi;
	dfpara[1] = lam0;
	double dx[14] = {0.0};
	int flag = GA_derivative_type1(0, x, dx, dfpara);

	double H = 0.0;
	H = V_Dot(costate, dx, 7);
	H = H + lam0*TmaxNU/Ispg0NU*u;
	return H;
}

// 二次同伦
double GA_Hamilton_type2(const double* x, double lam0, double epsi)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// 计算推力方向单位矢量和推力幅值大小
	double norm_lambdaV = V_Norm2(&costate[3], 3); // 速度协态的范数
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // 开关函数
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/2.0/epsi;
	// 计算状态变量的导数
	double dfpara[2];
	dfpara[0] = epsi;
	dfpara[1] = lam0;
	double dx[14] = {0.0};
	int flag = GA_derivative_type2(0, x, dx, dfpara);

	double H = 0.0;
	H = V_Dot(costate, dx, 7);
	H = H + lam0*TmaxNU/Ispg0NU*(u - epsi*u*(1.0 - u));
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
	// 获取MJD_MARS时刻火星位置
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
	// FILE *fid1 = fopen("temp1.txt", "w");
	// FILE *fid2 = fopen("temp2.txt", "w");
	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	int flag,NumPoint;
	double work[140]={0.0};
	flag = ode45(GA_derivative, x0, dfpara,  0.0, tm*TOF*86400/TUnit, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid1);
	// 计算一部分偏差
	// 引力辅助位置约束,7-9
	double rv_tm[6];
	rv02rvf(flag, rv_tm, rv_mars, MJD_MARS*86400/TUnit, (MJD0+tm*TOF)*86400/TUnit, muNU);
	for (int i=0;i<3;++i)
		fvec[7+i] = x0[i] - rv_tm[i];
	// normvinfin==normvinfout,10
	double vinfin[3], vinfout[3], normvinfin, normvinfout, unit_in[3], unit_out[3];
	V_Minus(vinfin, &x0[3], &rv_tm[3], 3);
	V_Add(vinfout, vinfin, dvg, 3);
	normvinfin = V_Norm2(vinfin, 3);
	normvinfout = V_Norm2(vinfout, 3);
	V_Divid(unit_in, vinfin, normvinfin, 3);
	V_Divid(unit_out, vinfout, normvinfout, 3);
	fvec[10] = normvinfin - normvinfout;
	// 互补松弛条件,11
	double theta = acos(V_Dot(unit_in, unit_out, 3));
	double rp = mpp/(normvinfin*normvinfout) * (1/sin(theta/2) - 1)/rmin;
	fvec[11] = kappa*(1 - rp);
	// fvec[11] = 1 - rp/rmin;
	// fvec[11] = rmin - rp;
	// fvec[11] = rp - rmin;

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

	// tm-时刻速度协态,12-14
	for (int i=0;i<3;++i)
		fvec[12+i] = x0[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	double H1 =  GA_Hamilton(x0, lam0, epsi);
	// 计算新的状态变量和协态变量
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // 速度
		x0[7+i] = x0[7+i] - chi[i]; // 位置协态
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // 速度协态
	}
	double H2 = GA_Hamilton(x0, lam0, epsi);
	// 静态条件偏差,15
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - kappa*c;


	// 第二阶段的积分
	// flag = ode45(GA_derivative, x0, dfpara,  tm*TOF*86400/TUnit,  tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	flag = ode45(GA_derivative, x0, dfpara, tm*TOF*86400/TUnit, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);

	for (int i=0;i<6;++i)
		fvec[i] = x0[i] - rvf[i];
	fvec[6] = x0[13];
	fvec[16] = enorm(13, x) - 1.0;

	if(outflag>0)
	{
		fvec[0]=x0[6];
	}
	// fclose(fid1);
	// fclose(fid2);
	return 0;
}

int GA_fvec_type1(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if(iflag==0)
	{
		return 0;
	}
	
	double x0[14], xf[14], rvf[6];
	V_Copy(x0, sfpara, 7);
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	double epsi = sfpara[14];
	epsi = 0.0;
	int outflag = NINT(sfpara[15]);
	// 获取MJD_MARS时刻火星位置
	double rv_mars[6];
	V_Copy(rv_mars, &sfpara[16], 6);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7);
	double chi[4], kappa, dvg[3], tm;
	V_Copy(chi, &x[8], 4);
	kappa = x[12];
	V_Copy(dvg, &x[13], 3);
	tm = x[16];
	
	bang_integrate(x0, 0.0, tm*TOF*86400/TUnit, STEP, xf, 14, epsi, lam0);

	int flag;
	// 计算一部分偏差
	// 引力辅助位置约束,7-9
	double rv_tm[6];
	rv02rvf(flag, rv_tm, rv_mars, MJD_MARS*86400/TUnit, (MJD0+tm*TOF)*86400/TUnit, muNU);
	for (int i=0;i<3;++i)
		fvec[7+i] = xf[i] - rv_tm[i];
	// normvinfin==normvinfout,10
	double vinfin[3], vinfout[3], normvinfin, normvinfout, unit_in[3], unit_out[3];
	V_Minus(vinfin, &xf[3], &rv_tm[3], 3);
	V_Add(vinfout, vinfin, dvg, 3);
	normvinfin = V_Norm2(vinfin, 3);
	normvinfout = V_Norm2(vinfout, 3);
	V_Divid(unit_in, vinfin, normvinfin, 3);
	V_Divid(unit_out, vinfout, normvinfout, 3);
	fvec[10] = normvinfin - normvinfout;
	// 互补松弛条件,11
	double theta = acos(V_Dot(unit_in, unit_out, 3));
	double rp = mpp/(normvinfin*normvinfout) * (1/sin(theta/2) - 1)/rmin;
	fvec[11] = kappa*(1 - rp);
	// fvec[11] = 1 - rp/rmin;
	// fvec[11] = rmin - rp;
	// fvec[11] = rp - rmin;

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

	// tm-时刻速度协态,12-14
	for (int i=0;i<3;++i)
		fvec[12+i] = xf[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	V_Copy(x0, xf, 14);
	double H1 =  GA_Hamilton_type1(x0, lam0, epsi);
	// 计算新的状态变量和协态变量
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // 速度
		x0[7+i] = x0[7+i] - chi[i]; // 位置协态
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // 速度协态
	}
	double H2 = GA_Hamilton_type1(x0, lam0, epsi);
	// 静态条件偏差,15
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - kappa*c;


	// 第二阶段的积分
	// flag = ode45(GA_derivative, x0, dfpara,  tm*TOF*86400/TUnit,  tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	// flag = ode45(GA_derivative, x0, dfpara, tm*TOF*86400/TUnit, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	bang_integrate(x0, tm*TOF*86400/TUnit, tf, STEP, xf, 14, epsi, lam0);

	for (int i=0;i<6;++i)
		fvec[i] = xf[i] - rvf[i];
	fvec[6] = xf[13];
	fvec[16] = enorm(13, x) - 1.0;

	if(outflag>0)
	{
		fvec[0]=xf[6];
	}
	// fclose(fid1);
	// fclose(fid2);
	return 0;
}

int GA_fvec_type2(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if(iflag==0)
	{
		return 0;
	}
	
	double x0[14], xf[14], rvf[6];
	V_Copy(x0, sfpara, 7);
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	double epsi = sfpara[14];
	int outflag = NINT(sfpara[15]);
	// 获取MJD_MARS时刻火星位置
	double rv_mars[6];
	V_Copy(rv_mars, &sfpara[16], 6);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7);
	double chi[4], kappa, dvg[3], tm;
	V_Copy(chi, &x[8], 4);
	kappa = x[12];
	V_Copy(dvg, &x[13], 3);
	tm = x[16];
	
	bi_integrate(x0, 0.0, tm*TOF*86400/TUnit, STEP, xf, 14, epsi, lam0);

	int flag;
	// 计算一部分偏差
	// 引力辅助位置约束,7-9
	double rv_tm[6];
	rv02rvf(flag, rv_tm, rv_mars, MJD_MARS*86400/TUnit, (MJD0+tm*TOF)*86400/TUnit, muNU);
	for (int i=0;i<3;++i)
		fvec[7+i] = xf[i] - rv_tm[i];
	// normvinfin==normvinfout,10
	double vinfin[3], vinfout[3], normvinfin, normvinfout, unit_in[3], unit_out[3];
	V_Minus(vinfin, &xf[3], &rv_tm[3], 3);
	V_Add(vinfout, vinfin, dvg, 3);
	normvinfin = V_Norm2(vinfin, 3);
	normvinfout = V_Norm2(vinfout, 3);
	V_Divid(unit_in, vinfin, normvinfin, 3);
	V_Divid(unit_out, vinfout, normvinfout, 3);
	fvec[10] = normvinfin - normvinfout;
	// 互补松弛条件,11
	double theta = acos(V_Dot(unit_in, unit_out, 3));
	double rp = mpp/(normvinfin*normvinfout) * (1/sin(theta/2) - 1)/rmin;
	fvec[11] = kappa*(1 - rp);
	// fvec[11] = 1 - rp/rmin;
	// fvec[11] = rmin - rp;
	// fvec[11] = rp - rmin;

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

	// tm-时刻速度协态,12-14
	for (int i=0;i<3;++i)
		fvec[12+i] = xf[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	V_Copy(x0, xf, 14);
	double H1 =  GA_Hamilton_type2(x0, lam0, epsi);
	// 计算新的状态变量和协态变量
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // 速度
		x0[7+i] = x0[7+i] - chi[i]; // 位置协态
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // 速度协态
	}
	double H2 = GA_Hamilton_type2(x0, lam0, epsi);
	// 静态条件偏差,15
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - kappa*c;


	// 第二阶段的积分
	// flag = ode45(GA_derivative, x0, dfpara,  tm*TOF*86400/TUnit,  tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	// flag = ode45(GA_derivative, x0, dfpara, tm*TOF*86400/TUnit, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	bi_integrate(x0, tm*TOF*86400/TUnit, tf, STEP, xf, 14, epsi, lam0);

	for (int i=0;i<6;++i)
		fvec[i] = xf[i] - rvf[i];
	fvec[6] = xf[13];
	fvec[16] = enorm(13, x) - 1.0;

	if(outflag>0)
	{
		fvec[0]=xf[6];
	}
	// fclose(fid1);
	// fclose(fid2);
	return 0;
}

// 将16个随机猜测值转换为17个打靶变量
// [0]:lam0
// [1-7]:7,初始时刻协态变量
// [8-11]:4,等式约束乘子
// [12]:不等式约束乘子
// [13-15]:3,速度增量
// [16]:引力辅助时间
void guess2lam(const double* guessValue, double* lamValue)
{
	double temp1, temp2;
	int i;
	// 用1个猜测值作为前8个打靶变量的norm，然后根据协态归一化条件计算出5个乘子的norm
	double norm1, norm2;
	norm1 = guessValue[0];
	norm2 = sqrt(1 - norm1*norm1);
	// 用7个猜测值，结合norm1，产生前8个打靶变量
	double alpha1[7];
	V_Copy(alpha1, &guessValue[1], 7);
	for (i=0;i<3;++i)
		alpha1[i]  = alpha1[i]*0.5*DPI;
	alpha1[3] = (alpha1[3] - 0.5)*DPI;
	alpha1[4] = alpha1[4]*D2PI;
	alpha1[5] = (alpha1[5] - 0.5)*DPI;
	alpha1[6] = alpha1[6]*D2PI;
	
	lamValue[0] = norm1*sin(alpha1[0]);
	temp1 = norm1*cos(alpha1[0])*cos(alpha1[1]);
	temp2 = temp1*sin(alpha1[2]);
	lamValue[1] = temp2*cos(alpha1[3])*cos(alpha1[4]);
	lamValue[2] = temp2*cos(alpha1[3])*sin(alpha1[4]);
	lamValue[3] = temp2*sin(alpha1[3]);
	temp2 = temp1*cos(alpha1[2]);
	lamValue[4] = temp2*cos(alpha1[5])*cos(alpha1[6]);
	lamValue[5] = temp2*cos(alpha1[5])*sin(alpha1[6]);
	lamValue[6] = temp2*sin(alpha1[5]);
	
	lamValue[7] = norm1*cos(alpha1[0])*sin(alpha1[1]);
	// 用4个猜测值，结合norm2，产生5个约束的乘子
	double alpha2[4];
	V_Copy(alpha2, &guessValue[8], 4);
	alpha2[0] = alpha2[0]*0.5*DPI;
	alpha2[1] = (alpha2[1] - 0.5)*DPI;
	alpha2[2] = (alpha2[2] - 0.5)*DPI;
	alpha2[3] = alpha2[3]*D2PI;
	temp1 = norm2*cos(alpha2[0]);
	temp2 = temp1*cos(alpha2[1]);
	lamValue[8] = temp2*cos(alpha2[2])*cos(alpha2[3]);
	lamValue[9] = temp2*cos(alpha2[2])*sin(alpha2[3]);
	lamValue[10] = temp2*sin(alpha2[2]);
	lamValue[11] = temp1*sin(alpha2[1]);
	lamValue[12] = norm2*sin(alpha2[0]);
	// 用1个猜测值，产生速度增量幅值
	double vAmplitude = sqrt(mpp/rmin)*guessValue[12];
	// 用2个猜测值，结合vAmplitude，产生3个速度分量
	double phi, delta;
	phi = (guessValue[13] - 0.5)*DPI;
	delta = guessValue[14]*D2PI;
	lamValue[13] = vAmplitude*cos(phi)*cos(delta);
	lamValue[14] = vAmplitude*cos(phi)*sin(delta);
	lamValue[15] = vAmplitude*sin(phi);
	// 用1个猜测值，产生引力辅助时刻
	lamValue[16] = (2.0 + guessValue[15])*365.25/2201;
}

// 对数同伦 随机猜测
int GA_FOP(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t)
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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	int num = 0;
	int j;
	// double amp, phi, delta;
	while (num<MaxGuessNum)
	{
		for (j=0; j<16; ++j)
			guessArray[j] = (double)rand()/RAND_MAX;
		guess2lam(guessArray, x);
		/*
		for (j=0;j<17;++j)
			printf("%16.6e%s%\n",x[j],",");
		printf("\n");
		*/

		/*
		// 第一个结果，剩余质量为16033.297
		// 不满足不等式约束
		x[0] = 2.506074e-1;
		x[1] = -4.032907e-1;
		x[2] = -4.962522e-1;
		x[3] = 1.749647e-2;
		x[4] =  5.821978e-1;
		x[5] = -4.035206e-1;
		x[6] = -4.535971e-2;
		x[7] = 9.193989e-2;
		x[8] = 7.709933e-2;
		x[9] = 4.817207e-2;
		x[10] = 3.540463e-2;
		x[11] = -7.853713e-2;
		x[12] = -3.187694e-16;
		x[13] = 3.467385e-1;
		x[14] = 1.273974e-1;
		x[15] = 2.647024e-3;
		x[16] = 4.122586e-1;
		

		// 第二个结果，剩余质量为15360.170
		x[0] = 1.847764e-1;
		x[1] = -3.645439e-1;
		x[2] = -5.972262e-1;
		x[3] = -3.226376e-2;
		x[4] = 4.865135e-1;
		x[5] = -4.495845e-1;
		x[6] = -6.959225e-3;
		x[7] = 1.321809e-1;
		x[8] = -2.069497e-2;
		x[9] = -7.312004e-2;
		x[10] = 2.202347e-2;
		x[11] = -1.120832e-1;
		x[12] = 8.996530e-3;
		x[13] = 9.283842e-2;
		x[14] = 5.918095e-2;
		x[15] = -9.730948e-3;
		x[16] = 3.627381e-1;

		
		// epsi = 0.995
		x[0] = 1.854903e-1;
		x[1] = -3.644848e-1;
		x[2] = -5.971348e-1;
		x[3] = -3.231720e-2;
		x[4] = 4.864155e-1;
		x[5] = -4.494969e-1;
		x[6] = -7.000981e-3;
		x[7] = 1.322510e-1;
		x[8] = -2.070964e-2;
		x[9] = -7.318025e-2;
		x[10] = 2.207289e-2;
		x[11] = -1.122068e-1;
		x[12] = 9.006350e-3;
		x[13] = 9.282229e-2;
		x[14] = 5.920827e-2;
		x[15] = -9.728563e-3;
		x[16] = 3.627526e-1;
		*/
		// epsi=1.2e-5时的结果,16021.554
		x[0] = 8.152167179245462e-001;
		x[1] = -2.143450414873604e-001;
		x[2] = -2.975501858963105e-001;
		x[3] = -3.214474440438555e-002;
		x[4] = 2.826663298308946e-001;
		x[5] = -2.109591915546727e-001;
		x[6] = -8.493627886200947e-002;
		x[7] = 1.781235703614778e-001;
		x[8] = 2.541862128708982e-002;
		x[9] = -5.919471028395332e-002;
		x[10] = 7.655083822488865e-002;
		x[11] = -1.616604224358071e-001;
		x[12] = 2.044634272557275e-002;
		x[13] = 6.413659577823913e-002;
		x[14] = 9.066606371583727e-002;
		x[15] = -1.837023674860771e-003;
		x[16] = 3.878463163809591e-001;
		/*
		// epsi=2.3437500e-8,16021.641737147094000,fail
		x[0] = 8.1541730965308412e-001;
		x[1] = -2.1421302795905156e-001;
		x[2] = -2.9733850261287192e-001;
		x[3] = -3.2099774206881676e-002;
		x[4] = 2.825155619269602e-001;
		x[5] = -2.1081110962928706e-001;
		x[6] = -8.4964463129163184e-002;
		x[7] = 1.7812588504346041e-001;
		x[8] = 2.5477784651429049e-002;
		x[9] = -5.9131632797361514e-002;
		x[10] = 7.6568560564841071e-002;
		x[11] = -1.6166523785486842e-001;
		x[12] = 2.0459482727909314e-002;
		x[13] = 6.4102281588356805e-002;
		x[14] = 9.0687059332444989e-002;
		x[15] = -1.8315469802102997e-003;
		x[16] = 3.8786694703116403e-001;
		*/

		info = hybrd1(GA_fvec, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
		// std::cout << enorm(n,fvec) << std::endl;
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
	std::cout << num << std::endl;
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.15fkg\n", Out[0]*MUnit);
	// printf("lamda0为:%.6e\n", Out[1]);
	printf("17个打靶变量值为:\n");
	for (int j=1; j<=17; j++)
		printf("%.15e,\n", Out[j]);
	return flag;
}

// 对数同伦 同伦过程
int GA_FOP_2(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t)
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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	// epsi=1.0时的结果
	x[0] = 1.847764e-1;
	x[1] = -3.645439e-1;
	x[2] = -5.972262e-1;
	x[3] = -3.226376e-2;
	x[4] = 4.865135e-1;
	x[5] = -4.495845e-1;
	x[6] = -6.959225e-3;
	x[7] = 1.321809e-1;
	x[8] = -2.069497e-2;
	x[9] = -7.312004e-2;
	x[10] = 2.202347e-2;
	x[11] = -1.120832e-1;
	x[12] = 8.996530e-3;
	x[13] = 9.283842e-2;
	x[14] = 5.918095e-2;
	x[15] = -9.730948e-3;
	x[16] = 3.627381e-1;
	
	// 把上一个结果作为初值
	if (Out[0]>1e-8)
		V_Copy(x, &Out[1], 17);

	printf("当前同伦参数为:%.15e\n", epsi);
	info = hybrd1(GA_fvec, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
	if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
	{
		sfpara[15]=1.0;
		int _j = GA_fvec(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
		if(fvec[0]>0.0) // 剩余质量为正，且满足打靶精度要求，停止
		{
			flag=1;
			Out[0]=fvec[0];
			V_Copy(&Out[1], x, 17);
		}
		sfpara[15]=0.0;
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.15fkg\n", Out[0]*MUnit);
		// printf("lamda0为:%.6e\n", Out[1]);
		
		printf("17个打靶变量值为:\n");
		for (int j=1; j<=17; j++)
			printf("%.20e,\n", Out[j]);
		
	}
	else
		printf("求解失败%d\n",flag);

	return flag;
}

// bang-bang 用对数同伦结果作为初值
int GA_FOP_type1(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t)
{
	double sfpara[22] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rv1, 6);
	sfpara[13] = tof;
	sfpara[14] = 0.0;
	sfpara[15] = 0.0;
	V_Copy(&sfpara[16], rv_middle, 6);

	
	int info, flag = 0;
	const int n = 17;
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	// 对数同伦 epsi=1.2e-5时的结果
	x[0] = 8.152167179245462e-001;
	x[1] = -2.143450414873604e-001;
	x[2] = -2.975501858963105e-001;
	x[3] = -3.214474440438555e-002;
	x[4] = 2.826663298308946e-001;
	x[5] = -2.109591915546727e-001;
	x[6] = -8.493627886200947e-002;
	x[7] = 1.781235703614778e-001;
	x[8] = 2.541862128708982e-002;
	x[9] = -5.919471028395332e-002;
	x[10] = 7.655083822488865e-002;
	x[11] = -1.616604224358071e-001;
	x[12] = 2.044634272557275e-002;
	x[13] = 6.413659577823913e-002;
	x[14] = 9.066606371583727e-002;
	x[15] = -1.837023674860771e-003;
	x[16] = 3.878463163809591e-001;

	printf("bang-bang控制计算结果：\n");
	info = hybrd1(GA_fvec_type1, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
	if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
	{
		sfpara[15]=1.0;
		int _j = GA_fvec_type1(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
		if(fvec[0]>Out[0]) // 剩余质量为正，且满足打靶精度要求，停止
		{
			flag=1;
			Out[0]=fvec[0];
			V_Copy(&Out[1], x, 17);
		}
		sfpara[14]=0.0;
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
		// printf("lamda0为:%.6e\n", Out[1]);
		
		printf("17个打靶变量值为:\n");
		for (int j=1; j<=17; j++)
			printf("%.6e,\n", Out[j]);
		
	}
	else
		printf("求解失败%d\n",flag);

	return flag;
}

// 二次同伦 随机猜测
int GA_FOP_type2(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t)
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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	int num = 0;
	int j;
	// double amp, phi, delta;
	while (num<MaxGuessNum)
	{
		for (j=0; j<16; ++j)
			guessArray[j] = (double)rand()/RAND_MAX;
		guess2lam(guessArray, x);
		for (j=0;j<17;++j)
			printf("%16.6e%s%\n",x[j],",");
		printf("\n");

		/*
		// 第一组计算结果，不等式约束乘子为零，剩余质量为16633.139
		x[0] = 5.059988e-1;
		x[1] = -3.510439e-1;
		x[2] = -4.324779e-1;
		x[3] = 1.426043e-2;
		x[4] = 5.198075e-1;
		x[5] = -3.620174e-1;
		x[6] = -5.811510e-2;
		x[7] = 1.044576e-1;
		x[8] = 8.770364e-2;
		x[9] = 3.765305e-2;
		x[10] = 4.969229e-2;
		x[11] = -7.980374e-2;
		x[12] = 9.987540e-16;
		x[13] = 3.488134e-1;
		x[14] = 1.081177e-1;
		x[15] = -9.460807e-3;
		x[16] = 4.089086e-1;
		
		// 第二组计算结果，剩余质量为15735.957，认为是正确结果
		x[0] = 6.137291e-1;
		x[1] = -2.792804e-1;
		x[2] = -4.608278e-1;
		x[3] = -5.363877e-2;
		x[4] = 3.636849e-1;
		x[5] = -3.348964e-1;
		x[6] = -5.564393e-2;
		x[7] = 1.768778e-1;
		x[8] = -7.358387e-3;
		x[9] = -1.036135e-1;
		x[10] = 6.240425e-2;
		x[11] = -1.905573e-1;
		x[12] = 1.729442e-2;
		x[13] = 7.823747e-2;
		x[14] = 7.901755e-2;
		x[15] = -5.190058e-3;
		x[16] = 3.751106e-1;
		
		// 第三组计算结果，剩余质量为15458.937
		x[0] = 2.214995e-1;
		x[1] = 2.132710e-1;
		x[2] = 5.669776e-1;
		x[3] = 2.738201e-2;
		x[4] = -4.769494e-1;
		x[5] = 5.622211e-1;
		x[6] = -3.827401e-2;
		x[7] = 1.393004e-1;
		x[8] = -1.067973e-1;
		x[9] = -6.298917e-3;
		x[10] = 3.868794e-2;
		x[11] = -7.713973e-2;
		x[12] = -3.375318e-14;
		x[13] = 9.186601e-2;
		x[14] = 3.508144e-1;
		x[15] = 1.159120e-2;
		x[16] = 4.233791e-1;
		
		// 第四组计算结果，剩余质量为16287.891
		x[0] = 4.280059e-1;
		x[1] = -3.916981e-1;
		x[2] = -5.170270e-1;
		x[3] = -5.082789e-2;
		x[4] = 4.775190e-1;
		x[5] = -3.586965e-1;
		x[6] = -2.796766e-2;
		x[7] = 1.137957e-1;
		x[8] = -5.637174e-2;
		x[9] = -7.778244e-2;
		x[10] = 4.638942e-2;
		x[11] = -1.080929e-1;
		x[12] = -7.160789e-15;
		x[13] = 1.256892e-1;
		x[14] = 2.341013e-1;
		x[15] = -1.450502e-2;
		x[16] = 3.667574e-1;
		
		// epsi=0.085，剩余质量15999.912
		x[0] = 7.919180e-1;
		x[1] = -2.276970e-1;
		x[2] = -3.201659e-1;
		x[3] = -4.036419e-2;
		x[4] = 2.969932e-1;
		x[5] = -2.262774e-1;
		x[6] = -8.546398e-2;
		x[7] = 1.785882e-1;
		x[8] = 1.700684e-2;
		x[9] = -6.772456e-2;
		x[10] = 7.835782e-2;
		x[11] = -1.640158e-1;
		x[12] = 1.905718e-2;
		x[13] = 6.840902e-2;
		x[14] = 8.779804e-2;
		x[15] = -2.476456e-3;
		x[16] = 3.850029e-1;
		
		// epsi = 0.068,16005.793
		x[0] = 7.965290e-1;
		x[1] = -2.254639e-1;
		x[2] = -3.157990e-1;
		x[3] = -3.905240e-2;
		x[4] = 2.946150e-1;
		x[5] = -2.232579e-1;
		x[6] = -8.554327e-2;
		x[7] = 1.785128e-1;
		x[8] = 1.834173e-2;
		x[9] = -6.597443e-2;
		x[10] = 7.806397e-2;
		x[11] = -1.627025e-1;
		x[12] = 1.928125e-2;
		x[13] = 6.798229e-2;
		x[14] = 8.809942e-2;
		x[15] = -2.379520e-3;
		x[16] = 3.854693e-1;
		
		// epsi = 0.051, 16011.347
		x[0] = 8.018689e-1;
		x[1] = -2.225724e-1;
		x[2] = -3.107233e-1;
		x[3] = -3.740513e-2;
		x[4] = 2.915161e-1;
		x[5] = -2.197941e-1;
		x[6] = -8.548479e-2;
		x[7] = 1.784613e-1;
		x[8] = 2.003096e-2;
		x[9] = -6.405784e-2;
		x[10] = 7.766326e-2;
		x[11] = -1.617180e-1;
		x[12] = 1.958372e-2;
		x[13] = 6.732724e-2;
		x[14] = 8.855402e-2;
		x[15] = -2.252143e-3;
		x[16] = 3.860175e-1;
		
		// epsi = 0.0485,16012.120
		x[0] = 8.026757e-1;
		x[1] = -2.221156e-1;
		x[2] = -3.099497e-1;
		x[3] = -3.712857e-2;
		x[4] = 2.910276e-1;
		x[5] = -2.192707e-1;
		x[6] = -8.546468e-2;
		x[7] = 1.784481e-1;
		x[8] = 2.031323e-2;
		x[9] = -6.376615e-2;
		x[10] = 7.759613e-2;
		x[11] = -1.616141e-1;
		x[12] = 1.963177e-2;
		x[13] = 6.719702e-2;
		x[14] = 8.864391e-2;
		x[15] = -2.230608e-3;
		x[16] = 3.861101e-1;
		
		// epsi = 0.0473,16012.486
		x[0] = 8.030649e-1;
		x[1] = -2.218931e-1;
		x[2] = -3.095759e-1;
		x[3] = -3.699205e-2;
		x[4] = 2.907899e-1;
		x[5] = -2.190183e-1;
		x[6] = -8.545367e-2;
		x[7] = 1.784411e-1;
		x[8] = 2.045250e-2;
		x[9] = -6.362509e-2;
		x[10] = 7.756301e-2;
		x[11] = -1.615687e-1;
		x[12] = 1.965519e-2;
		x[13] = 6.713074e-2;
		x[14] = 8.868959e-2;
		x[15] = -2.219905e-3;
		x[16] = 3.861560e-1;
		
		// epsi = 0.0459,16012.909
		x[0] = 8.035205e-1;
		x[1] = -2.216309e-1;
		x[2] = -3.091378e-1;
		x[3] = -3.682952e-2;
		x[4] = 2.905099e-1;
		x[5] = -2.187230e-1;
		x[6] = -8.543960e-2;
		x[7] = 1.784325e-1;
		x[8] = 2.061825e-2;
		x[9] = -6.345963e-2;
		x[10] = 7.752360e-2;
		x[11] = -1.615197e-1;
		x[12] = 1.968284e-2;
		x[13] = 6.705013e-2;
		x[14] = 8.874505e-2;
		x[15] = -2.207090e-3;
		x[16] = 3.862106e-1;
		
		// epsi = 0.0443,16013.388
		x[0] = 8.040432e-1;
		x[1] = -2.213276e-1;
		x[2] = -3.086345e-1;
		x[3] = -3.663931e-2;
		x[4] = 2.901863e-1;
		x[5] = -2.183844e-1;
		x[6] = -8.542181e-2;
		x[7] = 1.784218e-1;
		x[8] = 2.081217e-2;
		x[9] = -6.326932e-2;
		x[10] = 7.747748e-2;
		x[11] = -1.614689e-1;
		x[12] = 1.971484e-2;
		x[13] = 6.695351e-2;
		x[14] = 8.881141e-2;
		x[15] = -2.191978e-3;
		x[16] = 3.862746e-1;
		
		// epsi = 0.0436,16013.596
		x[0] = 8.042725965700290e-001;
		x[1] = -2.211937537176773e-001;
		x[2] = -3.084134241041890e-001;
		x[3] = -3.655454262671428e-002;
		x[4] = 2.900435398918890e-001;
		x[5] = -2.182359764224797e-001;
		x[6] = -8.541341705144076e-002;
		x[7] = 1.784169164509153e-001;
		x[8] = 2.089857298400147e-002;
		x[9] = -6.318563171857791e-002;
		x[10] = 7.745691426350598e-002;
		x[11] = -1.614484860170150e-001;
		x[12] = 1.972898956529453e-002;
		x[13] = 6.690966138611852e-002;
		x[14] = 8.884147982213525e-002;
		x[15] = -2.185199480239959e-003;
		x[16] = 3.863031369778016e-001;
		
		// epsi = 0.0401,16014.620
		x[0] = 8.054248501109432e-001;
		x[1] = -2.205131621605458e-001;
		x[2] = -3.073009940382859e-001;
		x[3] = -3.611564767208081e-002;
		x[4] = 2.893190301396931e-001;
		x[5] = -2.174912851870250e-001;
		x[6] = -8.536517986960314e-002;
		x[7] = 1.783898736382844e-001;
		x[8] = 2.134575909842813e-002;
		x[9] = -6.276318498478462e-002;
		x[10] = 7.735028955156803e-002;
		x[11] = -1.613641434423965e-001;
		x[12] = 1.980103548141748e-002;
		x[13] = 6.667519088863626e-002;
		x[14] = 8.900179100053358e-002;
		x[15] = -2.149613934468473e-003;
		x[16] = 3.864512557736997e-001;
		*/
		// epsi = 0.0302,16017.362
		x[0] = 8.087134641697199e-001;
		x[1] = -2.184904341370915e-001;
		x[2] = -3.041123392009458e-001;
		x[3] = -3.471704743361110e-002;
		x[4] = 2.871824860565886e-001;
		x[5] = -2.153857839013599e-001;
		x[6] = -8.515251769732488e-002;
		x[7] = 1.782852718859805e-001;
		x[8] = 2.276884626652113e-002;
		x[9] = -6.152814675477464e-002;
		x[10] = 7.700384717291900e-002;
		x[11] = -1.613042751996191e-001;
		x[12] = 2.001685034740052e-002;
		x[13] = 6.585171660602766e-002;
		x[14] = 8.955785178225488e-002;
		x[15] = -2.029132260449991e-003;
		x[16] = 3.869258034013191e-001;
		
		info = hybrd1(GA_fvec_type2, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
		// std::cout << enorm(n,fvec) << std::endl;
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			int _j = GA_fvec_type2(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
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
	std::cout << num << std::endl;
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	// printf("lamda0为:%.6e\n", Out[1]);
	printf("17个打靶变量值为:\n");
	for (int j=1; j<=17; j++)
		printf("%.6e,\n", Out[j]);
	return flag;
}
// 二次同伦 同伦过程
int GA_FOP_type2_process(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t)
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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa的维数至少是544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	
	// epsi=1.0时的结果
	x[0] = 6.137291e-1;
	x[1] = -2.792804e-1;
	x[2] = -4.608278e-1;
	x[3] = -5.363877e-2;
	x[4] = 3.636849e-1;
	x[5] = -3.348964e-1;
	x[6] = -5.564393e-2;
	x[7] = 1.768778e-1;
	x[8] = -7.358387e-3;
	x[9] = -1.036135e-1;
	x[10] = 6.240425e-2;
	x[11] = -1.905573e-1;
	x[12] = 1.729442e-2;
	x[13] = 7.823747e-2;
	x[14] = 7.901755e-2;
	x[15] = -5.190058e-3;
	x[16] = 3.751106e-1;
	
	// 把上一个结果作为初值
	if (Out[0]>1e-8)
		V_Copy(x, &Out[1], 17);

	printf("当前同伦参数为:%.15f\n", epsi);
	info = hybrd1(GA_fvec_type2, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()的输出标志
	if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
	{
		sfpara[15]=1.0;
		int _j = GA_fvec_type2(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
		if(fvec[0]>0.0) // 剩余质量为正，且满足打靶精度要求，停止
		{
			flag=1;
			Out[0]=fvec[0];
			V_Copy(&Out[1], x, 17);
		}
		sfpara[15]=0.0;
		printf("求解成功%d\n",flag);
		printf("剩余质量为:%.15fkg\n", Out[0]*MUnit);
		// printf("lamda0为:%.6e\n", Out[1]);
		
		printf("17个打靶变量值为:\n");
		for (int j=1; j<=17; j++)
			printf("%.15e,\n", Out[j]);
		
	}
	else
		printf("求解失败%d\n",flag);

	return flag;
}

int change_epsi(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t)
{
	double epsi = 1.0;
	int flag;
	flag = GA_FOP_2(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_middle, PSO_t);
	while (1)
	{
		if (epsi>0.05)
			epsi = epsi - 0.005;
		else if (epsi > 0.005)
			epsi = epsi - 0.001;
		else if (epsi > 0.0011)
			epsi = epsi - 0.0001;
		else if (epsi > 0.00045)
			epsi = epsi - 0.00005;
		else if (epsi > 0.00011)
			epsi = epsi - 0.00001;
		else if (epsi > 0.000012)
			epsi = epsi - 0.000001;
		else if (epsi > 1e-8)
			epsi = epsi*0.5;
		else
			break;
		flag = GA_FOP_2(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_middle, PSO_t);
	}
	return flag;
}

int change_epsi_type2(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t)
{
	double epsi = 1.0;
	int flag;
	flag = GA_FOP_type2_process(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_middle, PSO_t);
	while (1)
	{
		if (epsi>0.085)
			epsi = epsi - 0.005;
		else if (epsi > 0.068)
			epsi = epsi - 0.001;
		else if (epsi > 0.0485)
			epsi = epsi - 0.0005;
		else if (epsi > 0.045)
			epsi = epsi - 0.0001;
		else if (epsi > 0.0226)
			epsi = epsi - 0.00005;
		else if (epsi > 0.00003)
			epsi = epsi - 0.00001;
		else if (epsi > 0.000003)
			epsi = epsi - 0.000001;
		else
			break;
		flag = GA_FOP_type2_process(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_middle, PSO_t);
	}
	return flag;
}