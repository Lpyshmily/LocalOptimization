#include "ExternHeader.h"
#include "Constant.h"
#include "gravity_assist.h"

void FlybyDvA(double* vin, double* vpl, double* vout, double* deltavg, double* dv, double& delta, double& rp)
{
	//  引力甩摆子程序，采用模型，动力甩摆，在引力辅助后，飞出借力行星影响球后开启发动机。
	//	输入参数：
	//	vin[3]:甩摆前日心速度
	//	vpl[3]:行星绕太阳速度
	//	vout[3]:甩摆后日心速度
	//	输出参数
	//	deltavg[3]:引力辅助产生的速度增量
	//	dv[3]:甩摆前后总的速度增量
	//	delta:甩摆参数，rad
	//	rp:甩摆半径
    double vinf[3];
	V_Minus(vinf, vin, vpl, 3);
	V_Minus(dv, vout, vin, 3);

	if(V_Dot(vinf, vinf, 3) <= 0.0 || abs(V_Dot(vinf, vinf, 3)) < 1.0E-12)
	{
		deltavg[0] = 0.0;
		deltavg[1] = 0.0;
		deltavg[2] = 0.0;
		delta = 0.0;
		rp = 0.001;
		return;
	}

	// 计算三个方向的单位矢量
	double unit1[3], unit2[3], unit3[3];
	double normvinf = V_Norm2(vinf, 3);
	V_Divid(unit1, vinf, normvinf, 3);
	V_Cross(unit3, vin, vpl);
	double temp = V_Norm2(unit3, 3);
	V_Divid(unit3, unit3, temp, 3);
	V_Cross(unit2, unit3, unit1);

	double lam1 = (V_Dot(dv, unit1, 3)) / normvinf;
	double lam2 = (V_Dot(dv, unit2, 3)) / normvinf;
	double lam3 = (V_Dot(dv, unit3, 3)) / normvinf;
	// double lam1 = -(V_Dot(dv, unit1, 3)) / normvinf;
	// double lam2 = -(V_Dot(dv, unit2, 3)) / normvinf;
	// double lam3 = -(V_Dot(dv, unit3, 3)) / normvinf;
	double fai = atan2(lam2, lam3);
	while(fai < 0.0)
		fai = fai + D2PI;
	fai = fmod(fai, D2PI);
	
	temp = sqrt((1.0 + lam1)*(1.0 + lam1) + lam2*lam2 + lam3*lam3);
	double gamma = acos((1.0 + lam1) / temp);

	double deltam = 2.0 * asin(1.0 / (1.0 + normvinf*normvinf*rmin/mpp));
	if(gamma <= deltam)
	{
		delta = gamma;
		rp = mpp / (normvinf*normvinf) * (1.0 / sqrt(0.5 * (1.0 - (1.0 + lam1) / temp)) - 1.0);
	}
	else
	{
		delta = deltam;
		rp = rmin;
	}

	deltavg[0] = normvinf * (sin(delta) * (cos(fai) * unit3[0] + sin(fai) * unit2[0]) + (cos(delta) - 1.0) * unit1[0]);
	deltavg[1] = normvinf * (sin(delta) * (cos(fai) * unit3[1] + sin(fai) * unit2[1]) + (cos(delta) - 1.0) * unit1[1]);
    deltavg[2] = normvinf * (sin(delta) * (cos(fai) * unit3[2] + sin(fai) * unit2[2]) + (cos(delta) - 1.0) * unit1[2]);

	return;
}

void FlybyDvB(double* vin, double* vpl, double* vout, double& dvt, double* dv, double& delta, double& rp)
{
	//引力甩摆子程序，采用模型，动力甩摆，在公共近星点开启发动机。
	//	输入参数：
	//	vin[3]:甩摆前日心速度
	//	vpl[3]:行星绕太阳速度
	//	vout[3]:甩摆后日心速度
	//	输出参数
	//	dvt:需要施加的脉冲大小(负值表示减速)，km/s
	//	dv[3]:甩摆前后总的速度增量，km/s
	//	delta:甩摆参数，rad
	//	rp:甩摆半径，km
	
	double vinfin[3], vinfout[3];
	V_Minus(vinfin, vin, vpl, 3);
	V_Minus(vinfout, vout, vpl, 3);
	V_Minus(dv, vout, vin, 3);

	double normvinfin2 = V_Dot(vinfin, vinfin, 3);
	double normvinfin = sqrt(normvinfin2);
	double normvinfout2 = V_Dot(vinfout, vinfout, 3);
	double normvinfout = sqrt(normvinfout2);

	delta = V_Dot(vinfin, vinfout, 3)/(normvinfin*normvinfout);
	delta = acos(delta);
	// 牛顿迭代
	double x1, x2, y1, y2, temp;
	x1 = 1e-6;
	x2 = 0.0;
	int num = 0;
	while ((fabs(x2-x1)>1e-8) && (num<100))
	{
		y1 = acos(-mpp/(mpp+x1*normvinfout2)) + acos(-mpp/(mpp+x1*normvinfin2)) - delta - DPI;
		y2 = acos(-mpp/(mpp+x2*normvinfout2)) + acos(-mpp/(mpp+x2*normvinfin2)) - delta - DPI;
		temp = (x1*y2 - x2*y1)/(y2 - y1);
		x2 = x1;
		x1 = temp;
		num = num + 1;
	}
	rp = temp;
	if (num==100)
	{
		dvt = 1e6;
		return;
	}

	if(rp>=rmin)
		dvt = sqrt(2*mpp/rp+normvinfout2) - sqrt(2*mpp/rp+normvinfin2);
	else
		dvt = 1e6;        //若甩摆位置低于甩摆高度，则将脉冲大小设为很大的值，在PSO搜索时自动舍弃
	return;
}
// para 17维
// 0~5 rv0
// 6~11 rv1
// 12~17 rv_middle
// 18 mjd_middle
double GA_PSO_obj(const double* px, const double* para)
{
	double x, rv0[6], rv1[6], rv_middle[6], mjd_middle;
	x = *px;
	V_Copy(rv0, para, 6);
	V_Copy(rv1, &para[6], 6);
	V_Copy(rv_middle, &para[12], 6);
	mjd_middle = para[18];
	
	double dt1 = TOF*x; // MJD
	double dt2 = TOF*(1-x); // MJD
	double middle_t = MJD0 + dt1; // 引力辅助的时刻
	// 获得引力辅助时刻(middle_t)火星的位置和速度
	int flag;
	double rv_t[6] = {0.0};
	rv02rvf(flag, rv_t, rv_middle, mjd_middle*86400/TUnit, middle_t*86400/TUnit, muNU);
	// 分别求解两个lambert问题，获得vin和vout
	double v1[3] = {0.0}, v2[3] = {0.0};
	double vin[3] = {0.0}, vout[3] = {0.0};
	double dv0[3] = {0.0}, dvf[3] = {0.0};
	double unith[3] = {0.0, 0.0, 1.0};
	double a, e;

	lambert(v1, v2, a, e, rv0, rv_t, dt1*86400/TUnit, unith, flag, muNU);
	V_Minus(dv0, v1, &rv0[3], 3);
	V_Copy(vin, v2, 3);

	lambert(v1, v2, a, e, rv_t, rv1, dt2*86400/TUnit, unith, flag, muNU);
	V_Minus(dvf, &rv1[3], v2, 3);
	V_Copy(vout, v1, 3);

	// 引力辅助后开发动机
	
	double deltavg[3], dv[3];
	double delta, rp;
	FlybyDvA(vin, &rv_t[3], vout, deltavg, dv, delta, rp);
	// std::cout << rp << std::endl;
	double dvt[3];
	V_Minus(dvt, dv, deltavg, 3);

	double m_dv0 = V_Norm2(dv0, 3);
	double m_dvt = V_Norm2(dvt, 3);
	double m_dvf = V_Norm2(dvf, 3);
	
	return m_dv0 + m_dvt + m_dvf;
	

	// 公共近心点开启发动机
	/*
	double dv[3];
	double dvt, delta, rp;
	FlybyDvB(vin, &rv_t[3], vout, dvt, dv, delta, rp);
	double m_dv0 = V_Norm2(dv0, 3);
	
	double m_dvf = V_Norm2(dvf, 3);
	
	return m_dv0 + dvt + m_dvf;
	*/
}

double GA_PSO(const double* rv0, const double* rv1, double mjd_middle, const double* rv_middle, int psotime)
{
	double sfpara[19];
	V_Copy(sfpara, rv0, 6);
	V_Copy(&sfpara[6], rv1, 6);
	V_Copy(&sfpara[12], rv_middle, 6);
	sfpara[18] = mjd_middle;
	
	const int D = 1;
	double fbest, xbest[D];
	double wa[100];
	double fBestEver = 1e10, xBestEver[D];
	for(int i=0; i<psotime; i++)
	{
		cout<<i+1<<endl;
		PSO(GA_PSO_obj, xbest, fbest, sfpara, D, 20, wa);		
		if(fbest <= fBestEver)
		{
			fBestEver = fbest;
			for(int j=0; j<D; j++)
				xBestEver[j] = xbest[j];
		}
	}
	printf("引力辅助时刻为%.3f，速度增量为%.6f\n", *xBestEver, fBestEver*VUnit);
	return *xBestEver;
}