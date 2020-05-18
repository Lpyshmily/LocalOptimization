#include"OrbitFun.h"
//namespace OrbitFun{
//#define DPI 3.1415926535897932384626433832795
//#define D2PI 2.0*PI

// *             清华大学航天航空学院动力学与控制研究室博士生               *
// *                蒋方华(jiangfh04@mails.thu.edu.cn)                      *
// *                 最近修改: 2009.2.24                                    *


// 根据经典轨道根数求地心惯性直角坐标系下的位置和速度分量
//输入参数 coe[5],mu:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,可将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,可将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,可将其值取为omega+f,即纬度幅角
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14

//输出参数RV:
//   位置列矢量R分量X、Y、Z和速度列矢量V的分量VX、VY、VZ的6维向量
//   单位：米、米每秒.
//成功flag返回1，否则0
void coe2rv(int& flag, double* rv, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<=0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)
		return;
	if((coe[1]*cos(coe[5]))<-1.0)
	{
//		cout<<"不可能达到的双曲轨道."<<endl;		
		return;
	}

	double p=coe[0]*fabs(1.0-coe[1]*coe[1]);//半通径
	if(coe[1]==1.0)//如果是抛物线轨道,区别对待.
		p=2.0*coe[0];	

	double sini, cosi, sinO, cosO, sino, coso;
	sini=sin(coe[2]);
	cosi=cos(coe[2]);
	sinO=sin(coe[3]);
	cosO=cos(coe[3]);
	sino=sin(coe[4]);
	coso=cos(coe[4]);

	//轨道平面法向单位矢量,即角动量单位矢量
	double HVector[3]={sini*sinO, -sini*cosO, cosi};
      
	//偏心率单位矢量,或叫Laplace矢量
	double PVector[3]={cosO*coso-sinO*sino*cosi, sinO*coso+cosO*sino*cosi, sino*sini};

	//半通径方向单位矢量,PVector,QVector,HVector构成右手坐标系
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r=0.0;
	if((coe[1]*cos(coe[5]))+1.0<=0.0)
	{
//		cout<<"抛物或双曲轨道到达无穷远处."<<endl;
		r=1.0e308;
	}
	else
		r=p/(1.0+coe[1]*cos(coe[5]));	

	for(int i=0;i<3;i++)
	{
		rv[i]=r*(cos(coe[5])*PVector[i]+sin(coe[5])*QVector[i]);
		rv[3+i]=sqrt(mu/p)*(-sin(coe[5])*PVector[i]+(cos(coe[5])+coe[1])*QVector[i]);
	}
	flag=1;
	return;
}

// E2f 根据偏近点角和偏心率求真近点角
//输入参数E, e:
//   E:偏近点角,单位：弧度,对于双曲轨道,指H,有r=a(ecoshH-1)
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//输出参数f:
// 	真近点角（弧度）
//   epsilon:最小精度,默认为1.0e-12
double E2f(int& flag, double E, double e)
{
	if(e<0.0) {flag=0;return E;}

	double f=0.0;
    if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double E0=fmod(E, D2PI);
		if(E0>DPI)
			E0-=D2PI;
		if(E0<-DPI)
			E0+=D2PI;
		f=2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*E0));
		f=f+E-E0;
	}
	else if(e>1.0)//双曲线轨道
		f=2.0*atan(sqrt((e+1.0)/(e-1.0))*tanh(0.5*E));
	else// (abs(e-1.0)<epsilon)抛物线轨道
	{
		f=E;
//		cout<<"抛物线轨道没有定义偏近点角.在此将真近点角的值置为E."<<endl;
	}
	flag=1;
	return f;	
}


// E2M 根据偏近点角和偏心率求平近点角
//输入参数E, e:
//   E:偏近点角,单位：弧度,对于双曲轨道,指H
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   epsilon:最小精度,默认为1.0e-12
//输出参数M:
// 	平近点角(弧度),对于双曲轨道,指N
double E2M(int& flag, double E, double e)
{
	if(e<0.0){flag=0;return E;}
	double M=0.0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double E0=fmod(E, D2PI);    
		M=E0-e*sin(E0);
		M=M+E-E0;
	}
	else if(e>1.0)//双曲线轨道
		M=e*sinh(E)-E;   
	else//(abs(e-1.0)<epsilon)抛物线轨道
	{
		M=E;
//		cout<<"抛物线轨道没有定义偏近点角.在此将平近点角值置为E."<<endl;
	}
	flag=1;
	return M;
}

// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
//输入参数f0,t,a,e,mu,MaxIter,epsilon:
//   f0:初始真近点角,单位：弧度
//   t:轨道深化时间,单位：秒
//   a:轨道半长轴，单位：米。对于抛物线轨道，指近星距，即p/2
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   MaxIter:最大迭代次数
//   epsilon:迭代绝对误差,默认为1e-14;
//输出参数ft:
// 	最终真近点角(弧度)
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if(mu<=0.0||MaxIter<1||a<=0.0||e<0.0){flag=0;return f0;}

	double ft=0.0;
	if((e>=0.0&&e<1.0)||(e>1.0))//圆,椭圆,双曲轨道
	{
		double E=f2E(flag,f0,e);
		if(flag==0) return f0;
		double M=E2M(flag,E,e);
		if(flag==0) return f0;
		M+=sqrt(mu/(a*a*a))*dt;
		E=M2E(flag,M,e,MaxIter,epsilon);
		if(flag==0) return f0;
		ft=E2f(flag,E,e);
		if(flag==0) return f0;
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		if((f0<-DPI)||(f0>DPI))
		{
//			cout<<"对于抛物线轨道，初始真近点角应在-180至180度之间."<<endl;
			flag=0;return f0;
		}
		else if(f0>DPI||f0<-DPI)
			ft=f0;
		else
		{
			double B=0.75*sqrt(2.0*mu/(a*a*a))*dt+0.5*tan(0.5*f0)*((tan(0.5*f0))*(tan(0.5*f0))+3.0);
			double B1B=B+sqrt(1.0+B*B);
			double tanv=0.0;
			if(fabs(dt)<D2PI*sqrt((a*a*a)/mu)/1000.0)//推进时间为小量的情况
			{
				double A=pow(B1B, 2.0/3.0);
				tanv=2.0*A*B/(1.0+(1.0+A)*A);
			}
			else//不是小量的情况
			{
				double temp=pow(B1B, 1.0/3.0);
				tanv=temp-1.0/temp;
			}
			ft=2.0*atan(tanv);
		}
	}
	flag=1;
	return ft;
}


// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
//输入参数f0,ft,a,e,mu,AbsTol:
//   f0:初始真近点角,单位：弧度
//   ft:终端真近点角,单位：弧度
//   a:轨道半长轴，单位：米。对于抛物线轨道，指近星距，即p/2
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数t:
//   t:轨道推进时间,单位：秒
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu)
{
	flag=0;
	if(mu<=0.0||a<=0.0||e<0.0||ft<f0)
		return 0.0;

	double dt=0.0;
	if(e>=1.0)
	{
		double maxangle=DPI-acos(1.0/e);
		if((Min(f0,ft)<-maxangle)||(Max(f0,ft)>maxangle))
		{
//			cout<<"不可能达到的双曲或抛物线轨道"<<endl;			
			return 0.0;
		}
		else if(f0<-maxangle||f0>maxangle||ft<-maxangle||ft>maxangle)
		{
//			cout<<"所需时间很难计算准确,可能为无穷,在此置为1.0e100.";
			dt=1.0e308;			
			return dt;
		}		
	}

	double omega=sqrt(mu/(a*a*a));
	double delta=0.0;
	if((e>=0.0&&e<1.0)||(e>1.0))
	{
		double E=f2E(flag,f0,e);
		double M0=E2M(flag,E,e);
		E=f2E(flag,ft,e);
		double Mt=E2M(flag,E,e);
		if(flag==0) return 0.0;
		delta=Mt-M0;
	}
	else// if(fabs(e-1.0)<epsilon)
	{
		double B1=tan(0.5*f0)*((tan(0.5*f0))*(tan(0.5*f0))+3.0);
		double B2=tan(0.5*ft)*((tan(0.5*ft))*(tan(0.5*ft))+3.0);
		delta=sqrt(2.0)/3.0*(B2-B1);
	}	
	dt=delta/omega;
	flag=1;
	return dt;
}


// f2E 根据真近点角和偏心率求偏近点角
//输入参数f, e:
//   f:真近点角,单位：弧度
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   epsilon:最小精度,默认为1.0e-12
//输出参数E:
// 	偏近点角（弧度）,对于双曲轨道,指H,有r=a(ecoshH-1)
double f2E(int& flag, double f, double e)
{
	if(e<0.0){flag=0;return f;}

	double E=0.0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double f0=fmod(f, D2PI);
		if(f0>DPI)
			f0-=D2PI;
		if(f0<-DPI)
			f0+=D2PI;
		E=2.0*atan(sqrt((1.0-e)/(1.0+e))*tan(0.5*f0));
		E+=f-f0;
	}
	else if(e>1.0)//双曲线轨道
	{
		if(f>DPI-acos(1.0/e)||f<-DPI+acos(1.0/e))
		{
//			cout<<"不可能达到的双曲轨道."<<endl;
			flag=0;
			return f;
		}
		else
			E=2.0*atanh(sqrt((e-1.0)/(1.0+e))*tan(0.5*f));
	}
	else//if(abs(e-1.0)<epsilon)抛物线轨道
	{
		E=f;
//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为f."<<endl;
	}
	flag=1;
	return E;
}


// Inertia2LVLH 根据主星的绝对位置R0和速度V0及从星的绝对位置R1和速度V1求从星在主星轨道坐标系中的相对位置和速度
//注:R0,V0,R1,V1都向惯性系中投影,r,v向主星轨道坐标系投影
//输入参数R0[2], V0[2], R1[2], V1[2], epsilon:
//   R0[2], R1[2]:主星和从星相对地心的位置列矢量分量X、Y、Z,单位：米
//   V0[2], V1[2]:主星和从星相对地心的速度列矢量分量VX、VY、VZ,单位：米每秒.
//   epsilon:最小精度,默认为1.0e-12
//输出参数rv[5]:从星在主星轨道坐标系中的相对位置和相对速度列矢量分量,单位：米,米每秒
void Inertia2LVLH(int& flag, double* rv, const double* RV0, const double* RV1)
{
	int i;
	
	double R0[3]={RV0[0], RV0[1], RV0[2]};
	double V0[3]={RV0[3], RV0[4], RV0[5]};
	double R1[3]={RV1[0], RV1[1], RV1[2]};
	double V1[3]={RV1[3], RV1[4], RV1[5]};
	double radius=V_Norm2(R0,3);//距离
	if(radius<=0.0){flag=0;return;}
	double UnitR[3];
	for(i=0;i<3;i++) UnitR[i]=R0[i]/radius;//径向单位矢量    
	double VectorH[3];
	V_Cross(VectorH, R0, V0);
	double h=V_Norm2(VectorH,3);//角动量值
	if(h<=0.0){flag=0;return;}
	double UnitH[3];
	for(i=0;i<3;i++) UnitH[i]=VectorH[i]/h;//轨道面法向单位矢量
	
	double omega[3];
	for(i=0;i<3;i++) omega[i]=VectorH[i]/(radius*radius);
	double UnitF[3];
	V_Cross(UnitF, UnitH, UnitR);
	double R_R[3];
	for(i=0;i<3;i++) R_R[i]=R1[i]-R0[i];
	double R_V[3];
	V_Cross(R_V,omega,R_R);
	for(i=0;i<3;i++) R_V[i]=V1[i]-V0[i]-R_V[i];
	rv[0]=V_Dot(UnitR,R_R,3);
	rv[1]=V_Dot(UnitF,R_R,3);
	rv[2]=V_Dot(UnitH,R_R,3);
	rv[3]=V_Dot(UnitR,R_V,3);
	rv[4]=V_Dot(UnitF,R_V,3);
	rv[5]=V_Dot(UnitH,R_V,3);
	flag=1;
	return;
}

// LVLH2Inertia 根据主星的绝对位置R0和速度V0及从星相对主星轨道坐标系中的位置r和速度v求从星的绝对位置和速度
//注:RV0,RV1向惯性系投影;rv向主星轨道坐标系投影
//输入参数RV0[5],rv[5],epsilon:
//   RV0[5]:主星相对地心的位置和速度列矢量分量X、Y、Z,单位：米,VX、VY、VZ,单位：米每秒
//   rv[5]:从星相对主星轨道坐标系中的位置和速度列矢量分量X、Y、Z,单位：米,VX、VY、VZ,单位：米每秒
//   epsilon:最小精度,默认为1.0e-12
//输出参数RV1[5]:从星相对地心的位置和速度列矢量分量X、Y、Z,单位：米,VX、VY、VZ,单位：米每秒
void LVLH2Inertia(int& flag, double* RV1, const double* RV0, const double* rv, double epsilon)
{
	int i;
	
	double R0[3]={RV0[0], RV0[1], RV0[2]};
	double V0[3]={RV0[3], RV0[4], RV0[5]};
	double r[3]={rv[0], rv[1], rv[2]};
	double v[3]={rv[3], rv[4], rv[5]};
	double radius=V_Norm2(R0, 3);//距离
	if(radius<=0.0){flag=0;return;}
	double UnitR[3];
	for(i=0;i<3;i++) UnitR[i]=R0[i]/radius;//径向单位矢量    
	double VectorH[3];
	V_Cross(VectorH, R0, V0);
	double h=V_Norm2(VectorH,3);//角动量值
	if(h<=0.0){flag=0;return;}
	double UnitH[3];
	for(i=0;i<3;i++) UnitH[i]=VectorH[i]/h;//轨道面法向单位矢量
	
	double omega[3];
	for(i=0;i<3;i++) omega[i]=VectorH[i]/(radius*radius);
	double UnitF[3];
	V_Cross(UnitF, UnitH, UnitR);
	double R1[3];
	for(i=0;i<3;i++) R1[i]=R0[i]+r[0]*UnitR[i]+r[1]*UnitF[i]+r[2]*UnitH[i];
	double R_V[3];
	for(i=0;i<3;i++) R_V[i]=R1[i]-R0[i];
	double tempv[3];
	V_Cross(tempv, omega, R_V);
	for(i=0;i<3;i++)
	{
		RV1[i]=R1[i];
		RV1[3+i]=V0[i]+tempv[i]+v[0]*UnitR[i]+v[1]*UnitF[i]+v[2]*UnitH[i];
	}
	flag=1;
	return;
}


// M2E 根据平近点角和偏心率求偏近点角
//输入参数M, e, MaxIter, epsilon:
//   M:平近点角,单位：弧度,对于双曲轨道,指N
//   e:偏心率,0圆,(0,1)椭圆,1抛物线,(1,...)双曲线
//   MaxIter:最大迭代次数,默认为60
//   epsilon:迭代绝对误差,默认为1e-14;
//输出参数E:
// 	偏近点角（弧度），对双曲轨道，指H
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if(epsilon<=0.0||MaxIter<1||e<0.0){flag=0;return M;}

	//迭代方法见《Solar System Dynamics》Chapter2,Carl D.Murray and Stanley F.Dermott著
	double E=0.0, Minus=0.0, DeMinus=0.0, DeDeMinus=0.0, DeDeDeMinus=0.0, Delta1=0.0, Delta2=0.0, Delta3=0.0;
	int N=0;
	if(e>=0.0&&e<1.0)//圆和椭圆轨道
	{
		double RM=fmod(M, D2PI);
		if(RM<0.0)
			RM+=D2PI;
		double sinRM=sin(RM);
		E=RM+0.85*e*Sign(sinRM);
		N=0;   
		Delta3=1.0;
		while(fabs(Delta3)>=epsilon&&N<MaxIter)
		{
			Minus=E-e*sin(E)-RM;
			DeMinus=1.0-e*cos(E);
			DeDeMinus=e*sin(E);
			DeDeDeMinus=e*cos(E);
			Delta1=-Minus/DeMinus;
			Delta2=-Minus/(DeMinus+0.5*Delta1*DeDeMinus);
			Delta3=-Minus/(DeMinus+0.5*Delta2*DeDeMinus+1.0/6.0*Delta2*Delta2*DeDeDeMinus);
			E=E+Delta3;
			N=N+1;
		}    
		E=E+M-RM;
	}
	else if(e>1.0)//双曲线轨道
	{
		E=asinh(M/e);
		Delta3=1.0;
		N=0;
		while(fabs(Delta3)>=epsilon&&N<MaxIter)
		{
			Minus=e*sinh(E)-E-M;
			DeMinus=e*cosh(E)-1.0;
			DeDeMinus=e*sinh(E);
			DeDeDeMinus=e*cosh(E);
			Delta1=-Minus/DeMinus;
			Delta2=-Minus/(DeMinus+0.5*Delta1*DeDeMinus);
			Delta3=-Minus/(DeMinus+0.5*Delta2*DeDeMinus+1.0/6.0*Delta2*Delta2*DeDeDeMinus);
			E=E+Delta3;
			N=N+1;
		}  
	}
	else //(abs(e-1.0)<epsilon)抛物线轨道
	{
		E=M;
//		cout<<"抛物线轨道没有定义偏近点角.在此将其值置为M."<<endl;
	}
	if(((e>=0.0&&e<1.0)||(e>1.0))&&fabs(Delta3)>=5.0*epsilon&&N>=MaxIter)
	{
//		cout<<"迭代不收敛,请降低精度epsilon或增加迭代次数限制."<<endl;
		flag=0;
		return M;
	}
	flag=1;
	return E;
}



// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
//输入参数RV[5], mu, epsilon:
//   RV:位置列矢量分量X、Y、Z,单位：米,和速度列矢量分量VX、VY、VZ,单位：米每秒.
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数coe[6]:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距q,p=2q.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,将其值取为omega+f,即纬度幅角
void rv2coe(int& flag, double* coe, const double* RV, double mu)
{
	int i;
	flag=0;
	if(mu<=0.0)
		return;

	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R,3);//距离
	double velocity=V_Norm2(V,3);//速度
	if(radius<=0.0||velocity<=0.0)
		return;
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//径向单位矢量    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector,unitR,unitV);
	double h=radius*velocity*V_Norm2(hvector,3);//角动量值
	if(h<=0.0)
		return;
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector,3);//轨道面法向单位矢量
	//偏心率矢量
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];
	coe[1]=V_Norm2(evector,3);//偏心率
	double p=h*h/mu;
	if(coe[1]==1.0)
		coe[0]=0.5*p;//抛物线轨道的近星距
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//半长轴
	bool judge=(coe[1]>0.0);
	double unite[3]={0.0};
	if(judge)
		for(i=0;i<3;i++) unite[i]=evector[i]/coe[1];//偏心率单位矢量
	coe[2]=acos(unith[2]);//轨道倾角

	double unitN[3]={-unith[1], unith[0], 0.0};//节线矢量,未归一化

	double temp[3];

	if(V_Norm2(unitN,3)==0.0)
	{
		coe[3]=0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if(!judge)
		{
			coe[4]=0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5]=atan2(unitR[1]*unith[2],unitR[0]);//真近点角
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4]=atan2(unite[1]*unith[2],unite[0]); //近星点幅角       
			coe[5]=atan2(V_Dot(unith,temp,3), V_Dot(unite,unitR,3));
		}
    }
	else
	{
		V_Cross(temp, unitN, unitR);
		coe[3]=atan2(unith[0],-unith[1]);
		coe[5]=atan2(V_Dot(unith,temp,3), V_Dot(unitN,unitR,3));
		if(!judge)
		{
			coe[4]=0.0;
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4]=atan2(V_Dot(unith,temp,3), V_Dot(unite,unitN,3));
			coe[5]=coe[5]-coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3]=fmod(coe[3], D2PI);
	if(coe[3]<0.0)
		coe[3]+=D2PI;
	coe[4]=fmod(coe[4], D2PI);
	if(coe[4]<0.0)
		coe[4]+=D2PI;
	coe[5]=fmod(coe[5], D2PI);
	if(coe[1]>=1.0)
	{
		if(coe[5]>DPI-acos(1.0/coe[1]))
			coe[5]-=D2PI;
		else if(coe[5]<-DPI+acos(1.0/coe[1]))
			coe[5]+=D2PI;
	}
	flag=1;
	return;
}

void rv02rvf(int& flag, double* rv1, const double* rv0, double t0, double t1, double GM)
{
	double coe[6];
	rv2coe(flag,coe,rv0,GM);
	if(flag==0) return;
	coe[5]=f0dt2ft(flag,coe[5],t1-t0,coe[0],coe[1],GM);
	if(flag==0) return;
	coe[5]=fmod(coe[5],D2PI);
	coe2rv(flag,rv1,coe,GM);
}

void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
			  const double* rv0, const double* rv1, double t, double GM)
{
	double vt0[6],am,vt1[3],v0[3],v1[3],a,e,unith[3],dtemp0,dtemp1,dvmin=1.0e10,dv;
	int i,Nmax,n,flag0;
	flag=0;
	V_Cross(unith,rv0,&rv0[3]);
	dtemp0=V_Norm2(unith,3);
	if(dtemp0>= EPSILON)	for(i=0;i<3;i++) unith[i]/=dtemp0;
	else {unith[0]=0.0;unith[1]=0.0;unith[2]=1.0;}
	for(i=0;i<3;i++) vt0[i]=rv0[i]-rv1[i];
	am=(V_Norm2(rv0,3)+V_Norm2(rv1,3)+V_Norm2(vt0,3))/4.0;//最小椭圆转移轨道半长轴
	Nmax=(int)floor(t/(D2PI*sqrt(am*am*am/GM)));//最多可能的转移圈数
	for(n=0;n<=Nmax;n++)//圈数从0到最大可能数穷举，保留特征速度最小的结果
	{
		for(int j=0;j<2;j++)//对于多圈问题，有左右分枝两个解
		{
			lambert(v0, v1, a, e, rv0, rv1, t, unith, flag0, GM, 0, n, j);
			if(flag0!=1||e>=1.0) continue;//排除双曲和抛物轨道
			for(i=0;i<3;i++) 
			{
				vt0[i]=v0[i]-rv0[3+i];
				vt1[i]=rv1[3+i]-v1[i];
			}			
			dtemp0=V_Norm2(vt0,3);
			dtemp1=V_Norm2(vt1,3);
			dv=dtemp0+dtemp1;
			if(dv<dvmin)
			{
				flag=1;
				N=n;
				branch=j;
				dvmin=dv;
				V_Copy(dv0,vt0,3);
				V_Copy(dv1,vt1,3);
				Mdv0=dtemp0;
				Mdv1=dtemp1;
			}
			if(n==0) break;
		}
	}
}

// coe2ee 根据经典轨道根数求天球轨道根数，对轨道倾角180度奇异
//输入参数 coe[5],mu:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,可将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,可将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,可将其值取为omega+f,即纬度幅角
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数ee[6]:
// 	ee[0]:半通径p（米）.
// 	ee[1]:e*cos(omega+OMEGA),以f(不要与真近点角混淆)表示
// 	ee[2]:e*sin(omega+OMEGA),以g表示
//	ee[3]:tan(i/2)*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	ee(5:tan(i/2)*sin(OMEGA),以k表示
//	ee[5]:omega+OMEGA+f,以L表示
void coe2ee(int&flag, double* ee, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)	
		return;	
	if((coe[1]*cos(coe[5]))<-1.0)
	//		cout<<"不可能达到的双曲轨道."<<endl;
		return;
	
	ee[0]=coe[0]*fabs(1.0-coe[1]*coe[1]);//半通径p
	if(coe[1]==1.0)//如果是抛物线轨道,区别对待.
		ee[0]=2.0*coe[0];
	ee[1]=coe[1]*cos(coe[4]+coe[3]);//f
	ee[2]=coe[1]*sin(coe[4]+coe[3]);//g
	double temp=tan(coe[2]/2.0);
	ee[3]=temp*cos(coe[3]);//h
	ee[4]=temp*sin(coe[3]);//k
	ee[5]=coe[4]+coe[3]+coe[5];
	ee[5]=fmod(ee[5],D2PI);
	if(ee[5]<0.0)
		ee[5]+=D2PI;
	flag=1;
	return ;
}


// coe2mee 根据经典轨道根数求修正天球轨道根数
//输入参数 coe[5],mu:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,可将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,可将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,可将其值取为omega+f,即纬度幅角
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数[mee[6],orbtype]:
// 	mee[0]:半通径p（米）.
// 	mee[1]:e*cos(omega+I*OMEGA),以f(不要与真近点角混淆)表示
// 	mee[2]:e*sin(omega+I*OMEGA),以g表示
//	mee[3]:(tan(i/2))^I*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	mee(5:(tan(i/2))^I*sin(OMEGA),以k表示
//	mee[5]:omega+I*OMEGA+f,以L表示
//   orbtype:顺行轨道1,逆行轨道-1
void coe2mee(int&flag, double* mee, int& orbtype, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)	
		return;	
	if((coe[1]*cos(coe[5]))<-1.0)	
//		cout<<"不可能达到的双曲轨道."<<endl;
		return;

	orbtype=1;
	if(coe[2]>DPI/2.0)
		orbtype=-1;

	mee[0]=coe[0]*fabs(1.0-coe[1]*coe[1]);//半通径p
	if(coe[1]==1.0)//如果是抛物线轨道,区别对待.
		mee[0]=2.0*coe[0];
	mee[1]=coe[1]*cos(coe[4]+orbtype*coe[3]);//f
	mee[2]=coe[1]*sin(coe[4]+orbtype*coe[3]);//g
	double temp=tan(coe[2]/2.0);
	if(orbtype==-1)
		temp=1.0/tan(coe[2]/2.0);
	mee[3]=temp*cos(coe[3]);//h
	mee[4]=temp*sin(coe[3]);//k
	mee[5]=coe[4]+orbtype*coe[3]+coe[5];
	mee[5]=fmod(mee[5], D2PI);
	if(mee[5]<0.0)
		mee[5]+=D2PI;
	flag=1;
	return ;
}

// dee_dt 求天球轨道根数ee在沿径向、横向、法向加速度Ad作用下随时间的变化率
//输入参数 Ad, ee[5],mu:
//   Ad[0]:径向加速度大小，
//   Ad[1]:横向加速度大小，
//   Ad[0]:法向加速度大小，单位m/s^2
// 	ee[0]:半通径p（米）.
// 	ee[1]:e*cos(omega+OMEGA),以f(不要与真近点角混淆)表示
// 	ee[2]:e*sin(omega+OMEGA),以g表示
//	ee[3]:tan(i/2)*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	ee(5:tan(i/2)*sin(OMEGA),以k表示
//	ee[5]:omega+OMEGA+f,以L表示
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//输出参数dee[6]:分别是ee[6]每个分量对应的随时间变化率
void dee_dt(int&flag, double* dee, const double* Ad, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	double w=1.0+f*cos(L)+g*sin(L);
	double s2=1.0+h*h+k*k;
	double C1=sqrt(p/mu);
	dee[0]=C1*2.0*p/w*Ad[1];
	dee[1]=C1*( sin(L)*Ad[0]+((1.0+w)*cos(L)+f)/w*Ad[1]-(h*sin(L)-k*cos(L))*g/w*Ad[2]);
	dee[2]=C1*(-cos(L)*Ad[0]+((1.0+w)*sin(L)+g)/w*Ad[1]+(h*sin(L)-k*cos(L))*f/w*Ad[2]);
	dee[3]=C1*s2*cos(L)/2.0/w*Ad[2];
	dee[4]=C1*s2*sin(L)/2.0/w*Ad[2];
	dee[5]=sqrt(mu*p)*(w/p)*(w/p)+C1/w*(h*sin(L)-k*cos(L))*Ad[2];
	flag=1;
	return;
}


// dmee_dt 求修正天球轨道根数mee在沿径向、横向、法向加速度Ad作用下随时间的变化率
//输入参数 Ad, mee[5],orbtype,mu:
//   Ad[0]:径向加速度大小，
//   Ad[1]:横向加速度大小，
//   Ad[0]:法向加速度大小，单位m/s^2
// 	mee[0]:半通径p（米）.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),以f(不要与真近点角混淆)表示
// 	mee[2]:e*sin(omega+orbtype*OMEGA),以g表示
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	mee[4]:(tan(i/2))^orbtype*sin(OMEGA),以k表示
//	mee[5]:omega+orbtype*OMEGA+f,以L表示
//   orbtype:轨道标识,顺行轨道1,逆行轨道-1,默认为顺行轨道
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//输出参数dmee[6]:分别是mee[6]每个分量对应的随时间变化率
void dmee_dt(int& flag, double* dmee, const double* Ad, const double* mee, int orbtype, double mu)
{
	flag=0;
	if(mu<=0.0||abs(orbtype)!=1||mee[0]<=0.0||(mee[3]*mee[3]+mee[4]*mee[4])>1.0)
		return;

	double p=mee[0], f=mee[1], g=mee[2], h=mee[3], k=mee[4], L=mee[5];
	double w=1.0+f*cos(L)+g*sin(L);
	double s2=1.0+h*h+k*k;
	double C1=sqrt(p/mu);
	dmee[0]=C1*2.0*p/w*Ad[1];
	dmee[1]=C1*( sin(L)*Ad[0]+((1.0+w)*cos(L)+f)/w*Ad[1]-(orbtype*h*sin(L)-k*cos(L))*g/w*Ad[2]);
	dmee[2]=C1*(-cos(L)*Ad[0]+((1.0+w)*sin(L)+g)/w*Ad[1]+(orbtype*h*sin(L)-k*cos(L))*f/w*Ad[2]);
	dmee[3]=C1*s2*cos(L)/2.0/w*Ad[2]*orbtype;
	dmee[4]=C1*s2*sin(L)/2.0/w*Ad[2];
	dmee[5]=sqrt(mu*p)*(w/p)*(w/p)+C1/w*(orbtype*h*sin(L)-k*cos(L))*Ad[2];
	flag=1;
	return;
}

// ee2coe 根据天球根数求经典轨道根数
//输入参数 ee[5],mu,epsilon:
// 	ee[0]:半通径p（米）.
// 	ee[1]:e*cos(omega+OMEGA),以f(不要与真近点角混淆)表示
// 	ee[2]:e*sin(omega+OMEGA),以g表示
//	ee[3]:tan(i/2)*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	ee(5:tan(i/2)*sin(OMEGA),以k表示
//	ee[5]:omega+OMEGA+f,以L表示
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数coe[6]:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距q,p=2q.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,将其值取为omega+f,即纬度幅角
void ee2coe(int& flag, double* coe, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<=0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	coe[1]=sqrt(f*f+g*g);
	if(coe[1]==1.0)
		coe[0]=0.5*p;//抛物线轨道的近星距
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//半长轴
	double temp=sqrt(h*h+k*k);
	coe[2]=2.0*atan(temp);
	if(temp<=0.0)
	{
		coe[3]=0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5]=L;//真近点角
		}
		else
		{
			coe[4]=atan2(g,f); //近星点幅角       
			coe[5]=L-coe[4];
		}
	}
	else
	{
		coe[3]=atan2(k,h);
		coe[5]=L-coe[3];
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			coe[4]=atan2(g*h-f*k,f*h+g*k);
			coe[5]=coe[5]-coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3]=fmod(coe[3], D2PI);
	if(coe[3]<0.0)
		coe[3]+=D2PI;
	coe[4]=fmod(coe[4], D2PI);
	if(coe[4]<0.0)
		coe[4]+=D2PI;
	coe[5]=fmod(coe[5], D2PI);
	if(coe[5]<0.0)
		coe[5]+=D2PI;
	if(coe[1]>=1.0)
	{
		if(coe[5]>DPI-acos(1.0/coe[1]))
			coe[5]-=D2PI;
		else if(coe[5]<-DPI+acos(1.0/coe[1]))
			coe[5]+=D2PI;
	}
	flag=1;
	return;
}


// ee2rv 根据天球轨道根数求地心惯性直角坐标系下的位置和速度分量
//输入参数 ee[5],mu:
// 	ee[0]:半通径p（米）.
// 	ee[1]:e*cos(omega+OMEGA),以f(不要与真近点角混淆)表示
// 	ee[2]:e*sin(omega+OMEGA),以g表示
//	ee[3]:tan(i/2)*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	ee(5:tan(i/2)*sin(OMEGA),以k表示
//	ee[5]:omega+OMEGA+f,以L表示
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//输出参数RV:
//   位置列矢量R分量X、Y、Z和速度列矢量V的分量VX、VY、VZ组成的一个6维列向量,
//   单位：米、米每秒.
void ee2rv(int&flag, double* rv, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<=0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	double h_=sqrt(p/mu);
	double n=h_/(1.0+f*cos(L)+g*sin(L));
	double s2=1.0+h*h+k*k;
	double hh_kk=h*h-k*k;
	double hk2=2.0*h*k;

	double r=n*h_*mu;
	rv[0]=r/s2*(cos(L)+hh_kk*cos(L)+hk2*sin(L));
	rv[1]=r/s2*(sin(L)-hh_kk*sin(L)+hk2*cos(L));
	rv[2]=2.0*r/s2*(h*sin(L)-k*cos(L));
	rv[3]=-1.0/h_/s2*(sin(L)+hh_kk*sin(L)-hk2*cos(L)+g-hk2*f+hh_kk*g);
	rv[4]=-1.0/h_/s2*(-cos(L)+hh_kk*cos(L)+hk2*sin(L)-f+hk2*g+hh_kk*f);
	rv[5]=2.0/h_/s2*(h*cos(L)+k*sin(L)+f*h+g*k);
	flag=1;
	return;
}

// mee2coe 根据修正天球轨道根数求经典轨道根数
//输入参数 mee[5],orbtype,mu,epsilon:
// 	mee[0]:半通径p（米）.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),以f(不要与真近点角混淆)表示
// 	mee[2]:e*sin(omega+orbtype*OMEGA),以g表示
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),以k表示
//	mee[5]:omega+orbtype*OMEGA+f,以L表示
//   orbtype:轨道标识,顺行轨道1,逆行轨道-1,默认为顺行轨道
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数coe[6]:
// 	coe[0]半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
//                      其值取为近星距q,p=2q.
// 	coe[1]偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
// 	coe[2]轨道倾角i（弧度）:范围0<=i<=180度.
//	coe[3]升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,将其值取为0
//	coe[4]近拱点幅角omega（弧度）:当偏心率为0时没有意义,将其值置为0
//	coe[5]真近点角f（弧度）:当偏心率为0时没有意义,将其值取为omega+f,即纬度幅角
void mee2coe(int&flag, double* coe, const double* mee, int orbtype, double mu)
{
	flag=0;
	if(mu<=0.0||abs(orbtype)!=1||mee[0]<=0.0)
		return;

	double p=mee[0], f=mee[1], g=mee[2], h=mee[3], k=mee[4], L=mee[5];
	coe[1]=sqrt(f*f+g*g);
	if(coe[1]==1.0)
		coe[0]=0.5*p;//抛物线轨道的近星距
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//半长轴
	double temp=sqrt(h*h+k*k);
	if(temp>1.0)
	{
//		cout<<"第4与5个轨道根数平方和不能大于1."<<endl;
		return;
	}
	if(orbtype==1)
		coe[2]=2.0*atan(temp);
	else
		coe[2]=2.0*(DPI/2.0-atan(temp));
	if(temp<=0.0)
	{
		coe[3]=0.0;//升交点赤经
//		cout<<"轨道倾角接近0或180度,升交点赤经容易奇异.在此将其置为零."<<endl;
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;//近星点幅角
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;        
			coe[5]=L;//真近点角
		}
		else
		{
			coe[4]=atan2(g,f); //近星点幅角       
			coe[5]=L-coe[4];
		}
	}
	else
	{
		coe[3]=atan2(k,h);
		coe[5]=L-orbtype*coe[3];
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;
//			cout<<"偏心率接近0,近星点幅角容易奇异.在此将其置为零."<<endl;
		}
		else
		{
			coe[4]=atan2(g*h-orbtype*f*k,f*h+orbtype*g*k);
			coe[5]=coe[5]-coe[4];
		}
	}
	//转换到[0,2pi)中
	coe[3]=fmod(coe[3], D2PI);
	if(coe[3]<0.0)
		coe[3]+=D2PI;
	coe[4]=fmod(coe[4], D2PI);
	if(coe[4]<0.0)
		coe[4]+=D2PI;
	coe[5]=fmod(coe[5], D2PI);
	if(coe[5]<0.0)
		coe[5]+=D2PI;
	if(coe[1]>=1.0)
	{
		if(coe[5]>DPI-acos(1.0/coe[1]))
			coe[5]-=D2PI;
		else if(coe[5]<-DPI+acos(1.0/coe[1]))
			coe[5]+=D2PI;
	}
	flag=1;
	return;
}


// mee2rv 根据修正天球轨道根数求地心惯性直角坐标系下的位置和速度分量
//输入参数 mee[5],orbtype,mu:
// 	mee[0]:半通径p（米）.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),以f(不要与真近点角混淆)表示
// 	mee[2]:e*sin(omega+orbtype*OMEGA),以g表示
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),以k表示
//	mee[5]:omega+orbtype*OMEGA+f,以L表示
//   orbtype:轨道标识,顺行轨道1,逆行轨道-1,默认为顺行轨道
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//输出参数RV:
//   位置列矢量R分量X、Y、Z和速度列矢量V的分量VX、VY、VZ,
//   单位：米、米每秒.
void mee2rv(int&flag, double* rv, const double* mee, int orbtype, double mu)
{
	flag=0;
	if(mu<=0.0||abs(orbtype)!=1||mee[0]<=0.0)
		return;	

	double p=mee[0], f=mee[1], g=mee[2], h=mee[3], k=mee[4], L=mee[5];
	double temp=sqrt(h*h+k*k);
	if(temp>1.0)
	{
//		cout<<"第4与5个轨道根数平方和不能大于1."<<endl;
		return;
	}
	
	
	double h_=sqrt(p/mu);
	double n=h_/(1.0+f*cos(L)+g*sin(L));
	double s2=1.0+h*h+k*k;
	double hh_kk=h*h-k*k;
	double hk2=2.0*h*k;

	double r=n*h_*mu;
	rv[0]=r/s2*(cos(L)+hh_kk*cos(L)+hk2*orbtype*sin(L));
	rv[1]=r/s2*(orbtype*sin(L)-orbtype*hh_kk*sin(L)+hk2*cos(L));
	rv[2]=2.0*r/s2*(h*sin(L)-orbtype*k*cos(L));
	rv[3]=-1.0/h_/s2*(sin(L)+hh_kk*sin(L)-hk2*orbtype*cos(L)+g-hk2*orbtype*f+hh_kk*g);
	rv[4]=-1.0/h_/s2*(-orbtype*cos(L)+orbtype*hh_kk*cos(L)+hk2*sin(L)-f*orbtype+hk2*g+orbtype*hh_kk*f);
	rv[5]=2.0/h_/s2*(h*cos(L)+orbtype*k*sin(L)+f*h+orbtype*g*k);
	flag=1;
	return;
}

// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求天球轨道根数，对轨道倾角180度时奇异
//输入参数RV[5], mu, epsilon:
//   RV:位置列矢量分量X、Y、Z,单位：米,和速度列矢量分量VX、VY、VZ,单位：米每秒.
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数ee[6]:
// 	ee[0]:半通径p（米）.
// 	ee[1]:e*cos(omega+OMEGA),以f(不要与真近点角混淆)表示
// 	ee[2]:e*sin(omega+OMEGA),以g表示
//	ee[3]:tan(i/2)*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	ee(5:tan(i/2)*sin(OMEGA),以k表示
//	ee[5]:omega+OMEGA+f,以L表示
void rv2ee(int&flag, double* ee, const double* RV, double mu)
{
	flag=0;
	if(mu<=0.0)
		return;
	int i;
	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R, 3);//距离
	double velocity=V_Norm2(V, 3);//速度	
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//径向单位矢量    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h=radius*velocity*V_Norm2(hvector, 3);//角动量值
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector, 3);//轨道面法向单位矢量
	
	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//偏心率矢量	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];	
	//半能径方向矢量,模为e
	double qvector[3];
	V_Cross(qvector,unith,evector);
	//定义的一个矢量
	double unitA[3];
	for(i=0;i<3;i++) unitA[i]=h/mu*V[i]-qvector[i];
	//基本关系式
	// evector=e*[cos(omega)*cos(OMEGA)-sin(omega)*sin(OMEGA)*cos(i),
	//            cos(omega)*sin(OMEGA)+sin(omega)*cos(OMEGA)*cos(i),
	//            sin(omega)*sin(i)];
	// qvector=e*[-sin(omega)*cos(OMEGA)-cos(omega)*sin(OMEGA)*cos(i),
	//            -sin(omega)*sin(OMEGA)+cos(omega)*cos(OMEGA)*cos(i),
	//            cos(omega)*sin(i)];
	// unitR=[cos(OMEGA)*cos(omega+f)-cos(i)*sin(OMEGA)*sin(omega+f),
	//        sin(OMEGA)*cos(omega+f)+cos(i)*cos(OMEGA)*sin(omega+f),
	//        sin(i)*sin(omega+f)];
	// unitA=[-cos(OMEGA)*sin(omega+f)-cos(i)*sin(OMEGA)*cos(omega+f),
	//        -sin(OMEGA)*sin(omega+f)+cos(i)*cos(OMEGA)*cos(omega+f),
	//         sin(i)*cos(omega+f)];
	//推导得出
	//evector[0]+qvector[1]=(1+cos(i))*e*cos(omega+OMEGA);
	//evector[1]-qvector[0]=(1+cos(i))*e*sin(omega+OMEGA);
	//unith[0]=(1+cos(i))*tan(i/2)*sin(OMEGA);
	//unith[1]=-(1+cos(i))*tan(i/2)*cos(OMEGA);
	// unitR[0]+unitA[1]=(1+cos(i))*cos(omega+OMEGA+f);
	// unitR[1]-unitA[0]=(1+cos(i))*sin(omega+OMEGA+f);

	ee[0]=h*h/mu;//p

	if(unith[2]+1.0<=0.0)
	{
//		cout<<"轨道倾角接近180度，不适合用天球轨道根数描述. "<<endl;
		return;
	}
	double cosiadd1=1.0+unith[2];
	ee[1]=(evector[0]+qvector[1])/cosiadd1;//f
	ee[2]=(evector[1]-qvector[0])/cosiadd1;//g
	ee[3]=-unith[1]/cosiadd1;//h
	ee[4]=unith[0]/cosiadd1;//k
	ee[5]=atan2(unitR[1]-unitA[0],unitR[0]+unitA[1]);//L
	//转换到[0,2pi)中
	ee[5]=fmod(ee[5], D2PI);
	if(ee[5]<0.0)
		ee[5]+=D2PI;
	flag=1;
	return;
}


// rv2mee 根据地心惯性直角坐标系下的位置和速度分量求修正天球轨道根数
//输入参数RV[5], mu, epsilon:
//   RV:位置列矢量分量X、Y、Z,单位：米,和速度列矢量分量VX、VY、VZ,单位：米每秒.
//   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
//      默认为地心引力系数3.98600441800e+14
//   epsilon:最小精度,默认为1.0e-12
//输出参数[mee[6],orbtype]:
// 	mee[0]:半通径p（米）.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),以f(不要与真近点角混淆)表示
// 	mee[2]:e*sin(omega+orbtype*OMEGA),以g表示
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),以h(不要与角动量大小混淆)表示
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),以k表示
//	mee[5]:omega+orbtype*OMEGA+f,以L表示
//   orbtype:顺行轨道1,逆行轨道-1
void rv2mee(int&flag, double* mee, int& orbtype, const double* RV, double mu, double epsilon)
{
	flag=0;
	if(mu<=0.0)
		return;
	int i;

	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R, 3);//距离
	double velocity=V_Norm2(V, 3);//速度
	assert(radius>=epsilon&&velocity>=epsilon);
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//径向单位矢量    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//切向单位矢量
	double hvector[3];
	V_Cross(hvector,unitR,unitV);
	double h=radius*velocity*V_Norm2(hvector, 3);//角动量值
	assert(h>=epsilon);
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector, 3);//轨道面法向单位矢量
	
	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//偏心率矢量	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];	
	//半能径方向矢量,模为e
	double qvector[3];
	V_Cross(qvector,unith,evector);
	//定义的一个矢量
	double unitA[3];
	for(i=0;i<3;i++) unitA[i]=h/mu*V[i]-qvector[i];
	//基本关系式
	// evector=e*[cos(omega)*cos(OMEGA)-sin(omega)*sin(OMEGA)*cos(i),
	//            cos(omega)*sin(OMEGA)+sin(omega)*cos(OMEGA)*cos(i),
	//            sin(omega)*sin(i)];
	// qvector=e*[-sin(omega)*cos(OMEGA)-cos(omega)*sin(OMEGA)*cos(i),
	//            -sin(omega)*sin(OMEGA)+cos(omega)*cos(OMEGA)*cos(i),
	//            cos(omega)*sin(i)];
	// unitR=[cos(OMEGA)*cos(omega+f)-cos(i)*sin(OMEGA)*sin(omega+f),
	//        sin(OMEGA)*cos(omega+f)+cos(i)*cos(OMEGA)*sin(omega+f),
	//        sin(i)*sin(omega+f)];
	// unitA=[-cos(OMEGA)*sin(omega+f)-cos(i)*sin(OMEGA)*cos(omega+f),
	//        -sin(OMEGA)*sin(omega+f)+cos(i)*cos(OMEGA)*cos(omega+f),
	//         sin(i)*cos(omega+f)];
	//推导得出
	//evector[0]+orbtype*qvector[1]=(1+orbtype*cos(i))*e*cos(omega+OMEGA);
	//orbtype*evector[1]-qvector[0]=(1+orbtype*cos(i))*e*sin(omega+OMEGA);
	//unith[0]=(1+orbtype*cos(i))*(tan(i/2))^orbtype*sin(OMEGA);
	//unith[1]=-(1+orbtype*cos(i))*(tan(i/2))^orbtype*cos(OMEGA);
	// unitR[0]+orbtype*unitA[1]=(1+orbtype*cos(i))*cos(omega+OMEGA+f);
	// unitR[1]-orbtype*unitA[0]=(1+orbtype*cos(i))*sin(omega+OMEGA+f);

	mee[0]=h*h/mu;//p
	orbtype=1;
	if(unith[2]<0.0)
		orbtype=-1;
	double cosiadd1=1.0+orbtype*unith[2];
	mee[1]=(evector[0]+orbtype*qvector[1])/cosiadd1;//f
	mee[2]=(orbtype*evector[1]-qvector[0])/cosiadd1;//g
	mee[3]=-unith[1]/cosiadd1;//h
	mee[4]=unith[0]/cosiadd1;//k
	mee[5]=atan2(orbtype*unitR[1]-unitA[0], unitR[0]+orbtype*unitA[1]);//L
	//转换到[0,2pi)中
	mee[5]=fmod(mee[5], D2PI);
	if(mee[5]<0.0)
		mee[5]+=D2PI;
	flag=1;
	return;
}

//}//end namespace OrbitFun
