#include "integrate.h"
#include "ExternHeader.h"

// 计算开关函数及其一阶导数、二阶导数
// 返回值为开关函数的值， dSF:一阶导，ddSF:二阶导
double computeSF(const double* x, double& dSF, double& ddSF, double epsi, double lam0)
{
	double normlv = V_Norm2(&x[10], 3);
	double SF = 1.0 - (Ispg0NU*normlv/x[6] + x[13]) / lam0;
	double u;
	if (SF > epsi)
		u = 0.0;
	else if (SF < -epsi)
		u = 1.0;
	else
		u = 0.5*(1-SF/epsi);

	double C = Ispg0NU/lam0;
	double Dotlrlv = V_Dot(&x[7], &x[10], 3);
	dSF = C*Dotlrlv/x[6]/normlv;

	double normlr = V_Norm2(&x[7], 3);
	double radius = V_Norm2(x, 3);
	double C2 = muNU/(radius*radius*radius);
	double Dotrlv = V_Dot(x, &x[10], 3);
	ddSF = C/(normlv*x[6])*(-normlr*normlr + C2*normlv*normlv
		- 3.0*C2/(radius*radius)*Dotrlv*Dotrlv + (Dotlrlv/normlv)*(Dotlrlv/normlv)
		+ TmaxNU*u/Ispg0NU*Dotlrlv/x[6]);
	return SF;
}
// 二次同伦 满推
void GA_derivative_type2_1(double t, const double* x, double* dx, const double* dfpara)
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
	double u = 1.0;

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
}
// 中间推力
void GA_derivative_type2_2(double t, const double* x, double* dx, const double* dfpara)
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
}
// 无推力
void GA_derivative_type2_3(double t, const double* x, double* dx, const double* dfpara)
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

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3个位置分量的导数是3个速度分量
		dx[i+3] = -muNU/(r*r*r)*rv[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = 0.0;
	dx[13] = 0.0;
}

void bi_integrate(const double* x0, double t0, double tf, double h, double* xf, int dim, double epsi, double lam0)
{
	double SF0 = 1.0, SF1 = 1.0, dSF, ddSF;
	double dt = h, dt1 = 0.0, dt2 = 0.0;
	int i;
	double* xs = new double[dim];
	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = lam0;

	V_Copy(xs, x0, dim);
	while (t0 < tf)
	{
		dt = h;
		if (t0 + h > tf)
			dt = tf - t0;
		SF0 = computeSF(xs, dSF, ddSF, epsi, lam0);
		SF1 = SF0 + dSF*dt + 0.5*ddSF*dt*dt;
		if (SF0 <= -epsi)
		{
			if (SF1 <= -epsi)
				RK4(GA_derivative_type2_1, xs, t0, dt, xf, dim, dfpara); // 始终保持满推
			else
			{
				dt1 = 2.0*(SF0 + epsi)/(-dSF - sqrt(dSF*dSF - 2.0*ddSF*(SF0 + epsi)));
				RK4(GA_derivative_type2_1, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				if (SF1 <= epsi)
					RK4(GA_derivative_type2_2, xs, t0+dt1, dt - dt1, xf, dim, dfpara);
				else
				{
					dt2 = 2.0*(SF0 - epsi)/(-dSF - sqrt(dSF*dSF - 2.0*ddSF*(SF0 - epsi)));
					RK4(GA_derivative_type2_2, xs, t0+dt1, dt2 - dt1, xf, dim, dfpara);
					V_Copy(xs, xf, dim);
					RK4(GA_derivative_type2_3, xs, t0+dt2, dt - dt2, xf, dim, dfpara);
				}
			}
		}
		else if (SF0 > epsi)
		{
			if (SF1 > epsi)
				RK4(GA_derivative_type2_3, xs, t0, dt, xf, dim, dfpara);
			else
			{
				dt1 = 2.0*(SF0 - epsi)/(-dSF + sqrt(dSF*dSF - 2.0*ddSF*(SF0 - epsi)));
				RK4(GA_derivative_type2_3, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				if (SF1 >= -epsi)
					RK4(GA_derivative_type2_2, xs, t0 + dt1, dt - dt1, xf, dim, dfpara);
				else
				{
					dt2 = 2.0*(SF0 + epsi)/(-dSF + sqrt(dSF*dSF - 2.0*ddSF*(SF0 + epsi)));
					RK4(GA_derivative_type2_2, xs, t0 + dt1, dt2 - dt1, xf, dim, dfpara);
					V_Copy(xs, xf, dim);
					RK4(GA_derivative_type2_1, xs, t0 + dt2, dt - dt2, xf, dim, dfpara);
				}
			}
		}
		else
		{
			if (SF1 > epsi)
			{
				dt1 = 2.0*(SF0 - epsi)/(-dSF - sqrt(dSF*dSF - 2.0*ddSF*(SF0 - epsi)));
				RK4(GA_derivative_type2_2, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				RK4(GA_derivative_type2_3, xs, t0 + dt1, dt - dt1, xf, dim, dfpara);
			}
			else if (SF1 < -epsi)
			{
				dt1 = 2.0*(SF0 + epsi)/(-dSF + sqrt(dSF*dSF - 2.0*ddSF*(SF0 + epsi)));
				RK4(GA_derivative_type2_2, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				RK4(GA_derivative_type2_1, xs, t0 + dt1, dt - dt1, xf, dim, dfpara);
			}
			else
				RK4(GA_derivative_type2_2, xs, t0, dt, xf, dim, dfpara);
		}
		V_Copy(xs, xf, dim);
		t0 += dt;
	}
	delete[] xs;
}

double computeSF_type1(const double* x, double& dSF, double& ddSF, double epsi, double lam0)
{
	double normlv = V_Norm2(&x[10], 3);
	double SF = 1.0 - (Ispg0NU*normlv/x[6] + x[13]) / lam0;
	double u;
	if (SF > 0.0)
		u = 0.0;
	else
		u = 1.0;

	double C = Ispg0NU/lam0;
	double Dotlrlv = V_Dot(&x[7], &x[10], 3);
	dSF = C*Dotlrlv/x[6]/normlv;

	double normlr = V_Norm2(&x[7], 3);
	double radius = V_Norm2(x, 3);
	double C2 = muNU/(radius*radius*radius);
	double Dotrlv = V_Dot(x, &x[10], 3);
	ddSF = C/(normlv*x[6])*(-normlr*normlr + C2*normlv*normlv
		- 3.0*C2/(radius*radius)*Dotrlv*Dotrlv + (Dotlrlv/normlv)*(Dotlrlv/normlv)
		+ TmaxNU*u/Ispg0NU*Dotlrlv/x[6]);
	return SF;
}

// bang-bang 满推
void GA_derivative_type1_1(double t, const double* x, double* dx, const double* dfpara)
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
	double u = 1.0;

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
}

// 无推力
void GA_derivative_type1_2(double t, const double* x, double* dx, const double* dfpara)
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

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3个位置分量的导数是3个速度分量
		dx[i+3] = -muNU/(r*r*r)*rv[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = 0.0;
	dx[13] = 0.0;
}

// 在传入参数时，epsi任意
void bang_integrate(const double* x0, double t0, double tf, double h, double* xf, int dim, double epsi, double lam0)
{
	double SF0 = 1.0, SF1 = 1.0, dSF, ddSF;
	double dt = h, dt1 = 0.0, dt2 = 0.0;
	int i;
	double* xs = new double[dim];
	double dfpara[2] = {0.0};
	dfpara[0] = 0.0; // 强行设置同伦参数为0.0，在后面并没有用到，只是为了保持函数形式的一致
	dfpara[1] = lam0;

	V_Copy(xs, x0, dim);
	while (t0 < tf)
	{
		dt = h;
		if (t0 + h > tf)
			dt = tf - t0;
		SF0 = computeSF_type1(xs, dSF, ddSF, epsi, lam0);
		SF1 = SF0 + dSF*dt + 0.5*ddSF*dt*dt;
		if (SF0 <= 0.0)
		{
			if (SF1 <= 0.0)
				RK4(GA_derivative_type1_1, xs, t0, dt, xf, dim, dfpara); // 始终保持满推
			else
			{
				// dt1 = 2.0*(SF0 + epsi)/(-dSF - sqrt(dSF*dSF - 2.0*ddSF*(SF0 + epsi)));
				dt1 = 2.0*SF0/(-dSF - sqrt(dSF*dSF - 2.0*ddSF*SF0));
				RK4(GA_derivative_type1_1, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				RK4(GA_derivative_type1_2, xs, t0+dt1, dt - dt1, xf, dim, dfpara);
			}
		}
		else
		{
			if (SF1 > 0.0)
				RK4(GA_derivative_type1_2, xs, t0, dt, xf, dim, dfpara);
			else
			{
				// dt1 = 2.0*(SF0 - epsi)/(-dSF + sqrt(dSF*dSF - 2.0*ddSF*(SF0 - epsi)));
				dt = 2.0*SF0/(-dSF + sqrt(dSF*dSF - 2.0*ddSF*SF0));
				RK4(GA_derivative_type1_2, xs, t0, dt1, xf, dim, dfpara);
				V_Copy(xs, xf, dim);
				RK4(GA_derivative_type1_1, xs, t0+dt1, dt-dt1, xf, dim, dfpara);
			}
		}
		V_Copy(xs, xf, dim);
		t0 += dt;
	}
	delete[] xs;
}