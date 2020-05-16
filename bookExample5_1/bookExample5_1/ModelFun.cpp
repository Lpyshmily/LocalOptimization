#include "Constant.h"
#include "ExternHeader.h"
#include "ModelFun.h"


double L0Lt2dt(double L0, double Lt, const double* ee, double mus)
{
	int flag;
	double e = sqrt(ee[1]*ee[1] + ee[2]*ee[2]);
	double a = ee[0] / fabs(1.0 - e*e); // 保持a为正值
	double Oo = atan2(ee[2], ee[1]); // \omega + \Omega
	double dt = f0ft2dt(flag, L0-Oo, Lt-Oo, a, e, mus);
	return dt;
}

double L0dt2Lt(double L0, double dt, const double* ee, double mus)
{
	int flag;
	double e=sqrt(ee[1]*ee[1] + ee[2]*ee[2]);
	double a=ee[0] / fabs(1.0-e*e);
	double Oo=atan2(ee[2], ee[1]);
	double Lt=Oo + f0dt2ft(flag, L0-Oo, dt, a, e, mus);
	return Lt;
}

//M:6*3, dM/dx:6*3*6, dD6/dx:6, ddM/dx:6*3*6*6, ddD6/dx:6*6
//step:1,只算M, D6; 2, 算到dM, dD6;3, 算到ddM, ddD6
void MatMD(double * M, double& D6, double * dM, double* dD6, double* ddM, double* ddD6, const double* ee, int step)
{
	int i=0, j=0;

	double p=ee[0];
	double f=ee[1];
	double g=ee[2];
	double h=ee[3];
	double k=ee[4];
	double L=ee[5];

	double H=sqrt(p/muNU);
	double cosL=cos(L);
	double sinL=sin(L);
	double W=1.0+f*cosL+g*sinL;
	double S=1.0+h*h+k*k;
	double G=h*sinL-k*cosL;
	double HW=H/W;
	double HWG=HW*G;

	for(i=0;i<18;i++)
			M[i]=0.0;

	M[0*3+1]=2.0*p*HW;
	M[1*3+0]=H*sinL;
	M[1*3+1]=HW*((1.0+W)*cosL+f);
	M[1*3+2]=-HWG*g;
	M[2*3+0]=-H*cosL;
	M[2*3+1]=HW*((1.0+W)*sinL+g);
	M[2*3+2]=HWG*f;
	M[3*3+2]=0.5*HW*S*cosL;
	M[4*3+2]=0.5*HW*S*sinL;
	M[5*3+2]=HWG;
	D6=W*W/(H*H*H*muNU);
	if(step<2) 
		return;

	double cos2L=cos(2.0*L);
	double sin2L=sin(2.0*L);	
	double dWdL=-f*sinL+g*cosL;
	double dGdL=h*cosL+k*sinL;	
	double halfp=0.5/p;
	double cosLW=cosL/W;
	double sinLW=sinL/W;
	double dWdLW=dWdL/W;

	int idM;
		
	for(i=0;i<108;i++)
		dM[i]=0.0;			
	for(i=0;i<6;i++)
		dD6[i]=0.0;
	
	// M中共有18个元素，每个元素分别对六个轨道根数求导，M中的一行3个元素求导得18个元素，iDM=行数*18+列数*6作为M中每个元素求导的初始位置
	// M12
	idM=0*18+1*6;
	dM[idM]=3.0*HW;
	dM[idM+1]=-cosLW*M[0*3+1];
	dM[idM+2]=-sinLW*M[0*3+1];
	dM[idM+5]=-dWdLW*M[0*3+1];
	// M21
	idM=1*18+0*6;
	dM[idM]=halfp*M[1*3+0];
	dM[idM+5]=H*cosL;
	// M22
	idM=1*18+1*6;
	dM[idM]=halfp*M[1*3+1];
	dM[idM+1]=0.5*HW*(3.0+cos2L)-cosLW*M[1*3+1];
	dM[idM+2]=0.5*HW*sin2L-sinLW*M[1*3+1];
	dM[idM+5]=HW*(g*cos2L-f*sin2L-2.0*sinL)-dWdLW*M[1*3+1];
	idM=1*18+2*6;
	dM[idM]=halfp*M[1*3+2];
	dM[idM+1]=-cosLW*M[1*3+2];
	dM[idM+2]=-HWG-sinLW*M[1*3+2];
	dM[idM+3]=-HW*g*sinL;
	dM[idM+4]=HW*g*cosL;
	dM[idM+5]=-HW*g*dGdL-dWdLW*M[1*3+2];
	idM=2*18+0*6;
	dM[idM]=halfp*M[2*3+0];
	dM[idM+5]=H*sinL;
	idM=2*18+1*6;
	dM[idM]=halfp*M[2*3+1];
	dM[idM+1]=0.5*HW*sin2L-cosLW*M[2*3+1];
	dM[idM+2]=0.5*HW*(3.0-cos2L)-sinLW*M[2*3+1];
	dM[idM+5]=HW*(f*cos2L+g*sin2L+2.0*cosL)-dWdLW*M[2*3+1];
	idM=2*18+2*6;
	dM[idM]=halfp*M[2*3+2];
	dM[idM+1]=HWG-cosLW*M[2*3+2];
	dM[idM+2]=-sinLW*M[2*3+2];
	dM[idM+3]=HW*sinL*f;
	dM[idM+4]=-HW*cosL*f;
	dM[idM+5]=HW*f*dGdL-dWdLW*M[2*3+2];
	idM=3*18+2*6;
	dM[idM]=halfp*M[3*3+2];
	dM[idM+1]=-cosLW*M[3*3+2];
	dM[idM+2]=-sinLW*M[3*3+2];
	dM[idM+3]=HW*h*cosL;
	dM[idM+4]=HW*k*cosL;
	dM[idM+5]=-0.5*HW*S*sinL-dWdLW*M[3*3+2];
	idM=4*18+2*6;
	dM[idM]=halfp*M[4*3+2];
	dM[idM+1]=-cosLW*M[4*3+2];
	dM[idM+2]=-sinLW*M[4*3+2];
	dM[idM+3]=HW*h*sinL;
	dM[idM+4]=HW*k*sinL;
	dM[idM+5]=0.5*HW*S*cosL-dWdLW*M[4*3+2];
	idM=5*18+2*6;
	dM[idM+0]=halfp*HWG;
	dM[idM+1]=-cosLW*HWG;
	dM[idM+2]=-sinLW*HWG;
	dM[idM+3]=HW*sinL;
	dM[idM+4]=-HW*cosL;
	dM[idM+5]=HW*dGdL-dWdLW*HWG;

	double C=2.0*D6/W;
	dD6[0]=-1.5*D6/p;
	dD6[1]=C*cosL;
	dD6[2]=C*sinL;
	dD6[5]=C*dWdL;
	if(step<3)
		return;

	int q=0, z=0;
	int iddM;
	double ddWdL=1.0-W;
	double ddWdLW=ddWdL/W;
	double ddGdL=-G;
	for(i=0;i<648;i++)
		ddM[i]=0.0;			
	for(i=0;i<36;i++)
		ddD6[i]=0.0;
	//ddM[i*108+j*36+k*6+q]=ddM[i*108+j*36+q*6+k]
	
	idM=0*18+1*6+0;//第1行2列1纵指标索引
	iddM=0*108+1*36+0*6;
	C=-dM[idM];
	ddM[iddM]=1.5/p*HW;
	ddM[iddM+1]=C*cosLW;
	ddM[iddM+2]=C*sinLW;
	ddM[iddM+5]=C*dWdLW;
	iddM=0*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+2]*cosLW-dM[idM+1]*sinLW;
	ddM[iddM+5]=-dM[idM+5]*cosLW+M[0*3+1]*sinLW-dM[idM+1]*dWdLW;	
	iddM=0*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-M[0*3+1]*cosLW-dM[idM+5]*sinLW-dM[idM+2]*dWdLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-M[0*3+1]*ddWdLW;
	idM=1*18+0*6+0;//第2行1列1纵指标索引	
	iddM=1*108+0*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+5]=halfp*dM[idM+5];
	ddM[iddM+5*6+5]=-H*sinL;
	idM=1*18+1*6+0;//第2行2列1纵指标索引
	iddM=1*108+1*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+5]=halfp*dM[idM+5];	
	iddM=1*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW-HW*sin2L+M[1*3+1]*sinLW-dM[idM+5]*cosLW;
	iddM=1*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+2]*dWdLW+HW*cos2L-M[1*3+1]*cosLW-dM[idM+5]*sinLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-2.0*HW*(g*sin2L+f*cos2L+cosL)-M[1*3+1]*ddWdLW;
	idM=1*18+2*6+0;//第2行3列1纵指标索引
	iddM=1*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=1*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+2]*cosLW-dM[idM+1]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+5]*cosLW+M[1*3+2]*sinLW-dM[idM+1]*dWdLW;
	iddM=1*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-HW*sinL-dM[idM+3]*sinLW;
	ddM[iddM+4]=HW*cosL-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-HW*dGdL-M[1*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=1*108+2*36+3*6+5;
	C=HW*g;
	ddM[iddM]=-dM[idM+3]*dWdLW-C*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW-C*sinL;
	ddM[iddM+2*6]=-2.0*dM[idM+5]*dWdLW-C*ddGdL-M[1*3+2]*ddWdLW;
	idM=idM+0;//第3行1列1纵指标索引
	iddM=2*108+0*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+5]=halfp*dM[idM+5];
	ddM[iddM+5*6+5]=H*cosL;
	idM=2*18+1*6+0;//第3行2列1纵指标索引
	iddM=2*108+1*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=2*108+1*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+HW*cos2L+M[2*3+1]*sinLW-dM[idM+5]*cosLW;
	iddM=2*108+1*36+2*6+2;
	ddM[iddM]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+2]*dWdLW+HW*sin2L-M[2*3+1]*cosLW-dM[idM+5]*sinLW;
	ddM[iddM+3*6+3]=-2.0*dM[idM+5]*dWdLW-HW*2.0*(f*sin2L-g*cos2L+sinL)-M[2*3+1]*ddWdLW;
	idM=2*18+2*6+0;//第3行3列1纵指标索引
	iddM=2*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=2*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=HW*sinL-dM[idM+3]*cosLW;
	ddM[iddM+4]=-HW*cosL-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[2*3+2]*sinLW+HW*dGdL;
	iddM=2*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[2*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=2*108+2*36+3*6+5;
	ddM[iddM]=-dM[idM+3]*dWdLW+HW*f*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW+HW*f*sinL;
	ddM[iddM+2*6]=-2.0*dM[idM+5]*dWdLW+HW*ddGdL*f-M[2*3+2]*ddWdLW;
	idM=3*18+2*6+0;//第4行3列1纵指标索引
	iddM=3*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=3*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW-dM[idM+5]*cosLW+M[3*3+2]*sinLW;
	iddM=3*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[3*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=3*108+2*36+3*6+3;
	ddM[iddM]=HW*cosL;
	ddM[iddM+2]=-dM[idM+3]*dWdLW-HW*h*sinL;
	iddM=3*108+2*36+4*6+4;
	ddM[iddM]=HW*cosL;
	ddM[iddM+1]=-dM[idM+4]*dWdLW-HW*k*sinL;
	ddM[iddM+6+1]=-2.0*dM[idM+5]*dWdLW-0.5*HW*S*cosL-M[3*3+2]*ddWdLW;
	idM=4*18+2*6+0;//第5行3列1纵指标索引
	iddM=4*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=4*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[4*3+2]*sinLW-dM[idM+5]*cosLW;
	iddM=4*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[4*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=4*108+2*36+3*6+3;
	ddM[iddM]=HW*sinL;
	ddM[iddM+2]=-dM[idM+3]*dWdLW+HW*h*cosL;
	iddM=4*108+2*36+4*6+4;
	ddM[iddM]=HW*sinL;
	ddM[iddM+1]=-dM[idM+4]*dWdLW+HW*k*cosL;
	ddM[iddM+6+1]=-2.0*dM[idM+5]*dWdLW-0.5*HW*S*sinL-M[4*3+2]*ddWdLW;
	idM=5*18+2*6+0;//第6行3列1纵指标索引
	iddM=5*108+2*36+0*6;
	ddM[iddM]=-halfp*dM[idM];
	ddM[iddM+1]=halfp*dM[idM+1];
	ddM[iddM+2]=halfp*dM[idM+2];
	ddM[iddM+3]=halfp*dM[idM+3];
	ddM[iddM+4]=halfp*dM[idM+4];
	ddM[iddM+5]=halfp*dM[idM+5];
	iddM=5*108+2*36+1*6;
	ddM[iddM+1]=-2.0*dM[idM+1]*cosLW;
	ddM[iddM+2]=-dM[idM+1]*sinLW-dM[idM+2]*cosLW;
	ddM[iddM+3]=-dM[idM+3]*cosLW;
	ddM[iddM+4]=-dM[idM+4]*cosLW;
	ddM[iddM+5]=-dM[idM+1]*dWdLW+M[5*3+2]*sinLW-dM[idM+5]*cosLW;
	iddM=5*108+2*36+2*6;
	ddM[iddM+2]=-2.0*dM[idM+2]*sinLW;
	ddM[iddM+3]=-dM[idM+3]*sinLW;
	ddM[iddM+4]=-dM[idM+4]*sinLW;
	ddM[iddM+5]=-dM[idM+2]*dWdLW-M[5*3+2]*cosLW-dM[idM+5]*sinLW;
	iddM=5*108+2*36+3*6+5;
	ddM[iddM]=-dM[idM+3]*dWdLW+HW*cosL;
	ddM[iddM+6]=-dM[idM+4]*dWdLW+HW*sinL;
	ddM[iddM+12]=-2.0*dM[idM+5]*dWdLW+HW*ddGdL-M[5*3+2]*ddWdLW;
	for(i=0;i<6;i++)
	{
		for(j=0;j<3;j++)
			for(q=0;q<6;q++)
				for(z=0;z<q;z++)
					ddM[i*108+j*36+q*6+z]=ddM[i*108+j*36+z*6+q];
	}

	ddD6[0*6+0]=-2.5/p*dD6[0];
	C=-1.5/p;
	ddD6[0*6+1]=C*dD6[1];
	ddD6[0*6+2]=C*dD6[2];
	ddD6[0*6+5]=C*dD6[5];
	ddD6[1*6+0]=2.0*dD6[0]*cosLW;
	ddD6[1*6+1]=dD6[1]*cosLW;
	ddD6[1*6+2]=-dD6[1]*sinLW+2.0*dD6[2]*cosLW;
	ddD6[1*6+5]=-dD6[1]*dWdLW-2.0*D6*sinLW+2.0*dD6[5]*cosLW;
	ddD6[2*6+0]=2.0*dD6[0]*sinLW;
	ddD6[2*6+1]=-dD6[2]*cosLW+2.0*dD6[1]*sinLW;
	ddD6[2*6+2]=-dD6[2]*sinLW+2.0*dD6[2]*sinLW;
	ddD6[2*6+5]=-dD6[2]*dWdLW+2.0*D6*cosLW+2.0*dD6[5]*sinLW;
	ddD6[5*6+0]=2.0*dD6[0]*dWdLW;
	ddD6[5*6+1]=2.0*dD6[1]*dWdLW-dD6[5]*cosLW-2.0*D6*sinLW;
	ddD6[5*6+2]=2.0*dD6[2]*dWdLW-dD6[5]*sinLW+2.0*D6*cosLW;
	ddD6[5*6+5]=dD6[5]*dWdLW+2.0*D6*ddWdLW;	
}

//燃料最优 计算14个状态变量和协态变量的导数，即摄动方程
int DynEquFuelOpt(double t, const double* x, double* dx, const double* dfpara)
{
	double epsi, lam0,  m;
	double v_lrv[6]={0.0}, alpha[3]={0.0};
	double norma, u, rou, tempd, lm;
	int i, j, k;
	
	epsi=dfpara[0];
	lam0=dfpara[1];
		
	m=x[6];
	//if(m<0.2)
	//	return -1;
	lm=x[13];
	V_Copy(v_lrv, &x[7], 6);
	
	double M[18]={0.0};
	double dM[108]={0.0};
	double D6=0.0;
	double dD6[6]={0.0};
	double ddM[1]={0.0};
	double ddD6[1]={0.0};
	MatMD(M, D6, dM, dD6,  ddM, ddD6, x, 2);
	for(j=0;j<3;j++)
	{
		alpha[j]=0.0;
		for(i=0;i<6;i++)
			alpha[j]-=M[i*3+j]*v_lrv[i]; // 由极大值原理确定推力方向
	}
	norma=enorm(3, alpha);
	for(j=0;j<3;j++)
		alpha[j]/=norma;
	rou=1.0-(Ispg0NU*norma/m+lm)/lam0; // 开关函数
	// u=2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // 对数型同伦指标下的推力
	if (rou>0)
		u = 0;
	else
		u = 1;

	for(i=0;i<14;i++)
		dx[i]=0.0;
		
	 tempd=TmaxNU*u/m;
	 for(i=0;i<6;i++)
	 {
		 dx[i]=0.0;
		 for(j=0;j<3;j++)
			 dx[i]+=M[i*3+j]*alpha[j];
		 dx[i]*=tempd;
	 }
	 dx[5]+=D6;
	 dx[6]=-TmaxNU*u/Ispg0NU; // 质量的导数
	 for(i=0;i<6;i++)
	 {
		 dx[7+i]=0.0;
		 for(j=0;j<6;j++)
		 {
			 for(k=0;k<3;k++)
				 dx[7+i]-=v_lrv[j]*tempd*dM[j*18+k*6+i]*alpha[k];
		 }
		 dx[7+i]-=v_lrv[5]*dD6[i];
	 }
	 dx[13]=-TmaxNU*u/(m*m)*norma;	
	return 1;
}

//利用DynEquFuelOpt()得到14个积分变量随时间变化的导数，通过ode45()进行给定协态初值情况下的积分，得到估计的终端状态变量
//两端时间和状态固定的燃料最优问题打靶方程
//sfpara
//[0-6],初始春分点轨根和质量
//[7-12],末端春分点轨根
//[13-14],转移时间和epsilon
//[15],输出信息标识，0表示不输出，1表示输出
//x:[0],lam0;[1-7]初始协态
//n为方程维数，fvec为残差，iflag为0时表示上方调用函数告知不必算残差
//计算成功时，返回0或大于0的值，不成功时返回小于0的值
int ShootFunFuelOpt(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	double x0[14]={0.0}; // 初始状态变量和协态变量，各7个
	double AbsTol[14]={0.0};
	double dfpara[2]={0.0};
	double eef[6]={0.0}; // 末端状态变量，6个，不含质量
	double work[140]={0.0};

	double lam0=1.0, RelTol=1E-12, tf=0.0, epsi=1.0; // tf-终端时刻，epsi-同伦参数
	int i;
		
	V_Copy(x0, &sfpara[0], 7); // 赋予x0 7个状态变量初值
	V_Copy(eef, &sfpara[7], 6);	
	tf=sfpara[13];
	epsi=sfpara[14];
	int outflag=NINT(sfpara[15]);

	lam0=x[0];	
	V_Copy(&x0[7], &x[1], 7); // 赋予x0 7个协态变量初值
	
	dfpara[0]=epsi;
	dfpara[1]=lam0;
	
	if(iflag==0)
	{
		return 0;
	}
	for(i=0;i<14;i++)
		AbsTol[i]=1e-12;
	RelTol=1e-12;
	// FILE *fid=NULL;//fopen("temp.txt","w");//如果设定有效文件路径，最后需要关闭文件
	FILE *fid = fopen("temp.txt", "w");
	int flag,NumPoint;
	flag=ode45(DynEquFuelOpt, x0, dfpara,  0.0, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
	for(i=0;i<6;i++)
		fvec[i]=x0[i]-eef[i];
	fvec[6]=x0[13];
	fvec[7]=enorm(8, x)-1.0;
	if(outflag>0)
	{
		fvec[0]=x0[6];
	}
	fclose(fid);
	return 0;
}

//反复随机猜初值，通过ShootFunFuelOpt()计算给定初值条件下的积分结果，再通过hybrid()进行打靶方程求解
//求两端固定的燃料最优问题
//Out为输出结果，分别为剩余质量，8个乘子初值（第一个为lamd0)，共9维
//ee0为出发的春分点轨道根数，ee1为到达的春分点轨道根数，m0为出发质量，ft为飞行时间，
//epsi为同伦参数的值（直接取一个较小的值，如0.001，求得接近邦邦控制的结果），MaxGuessNum为最大猜测次数
//计算成功时，返回1，否则返回0
int Fixed2PFuelOpt(double* Out, const double* ee0, const double* ee1, double m0, double ft, double epsi, int MaxGuessNum)
{

	int j, n, info, flag=0;
	double xtol=1.0e-8;
	double x[8]={0.0}; // 8个协态变量初值，需要打靶求解
	double fvec[8]={0.0}, wa[200]={0.0};
	double sfpara[16]={0.0}; // 0~6-7个初始状态值，7~13-7个终端状态值，14-同伦参数，15-ShootFunFuelOpt()的输出信息标志
	V_Copy(sfpara, ee0, 6);
	sfpara[6]=m0;
	V_Copy(&sfpara[7], ee1, 6);
	sfpara[13]=ft;
	sfpara[14]=epsi;
	sfpara[15]=0.0;
	n = 8; // 打靶变量个数为8个
	Out[0]=0.0;
	int num=0;
	while(num<=MaxGuessNum)
	{
		for(j=0;j<8;j++)
				x[j]=(double)rand()/RAND_MAX-0.5; // 随机给出打靶变量初值
		x[7]+=0.5; // 质量协态导数为负，且末端为0，因此初值必为正
		x[0]+=0.5; // lambda0，归一化乘子，人为取正
		info = hybrd1(ShootFunFuelOpt, n, x, fvec, sfpara, wa, xtol, 20, 500); // info-hybrd1()的输出标志
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			j=ShootFunFuelOpt(n, x, fvec, 1, sfpara); // 利用同伦计算得到的协态初值，进行最后一次的积分求解（直接用这个较小的同伦参数，求解近似邦邦控制的结果）
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

void OutputAmplitude(const double* Out, const double* ee0, double m0, double epsi, const double* TimeNode, int TimeNodeNum)
{
	double x[14] = {0.0};
	V_Copy(x, ee0, 6);
	x[6] = m0;
	V_Copy(&x[7], &Out[2], 7);
	double original_x[14] = {0.0};
	V_Copy(original_x, x, 14);

	double dfpara[2] = {0.0};
	dfpara[0] = epsi;
	dfpara[1] = Out[1];

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	FILE *fid = NULL;
	FILE *ufid = fopen("amplitude.txt", "w");
	int flag,NumPoint;
	double work[140]={0.0};

	double M[18]={0.0};
	double dM[108]={0.0};
	double D6=0.0;
	double dD6[6]={0.0};
	double ddM[1]={0.0};
	double ddD6[1]={0.0};
	double alpha[3] = {0.0}; // 推力方向单位矢量
	double norma, rou, u;

	for (int i=0;i<TimeNodeNum;++i)
	{
		V_Copy(x, original_x, 14);
		
		// 积分到第i个时间节点
		flag = ode45(DynEquFuelOpt, x, dfpara, 0.0, TimeNode[i], 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid);
		// 计算第i个时间节点的推力大小
		MatMD(M, D6, dM, dD6,  ddM, ddD6, x, 2); // 计算M矩阵
		for(int j=0;j<3;j++)
		{
			alpha[j]=0.0;
			for(int k=0;k<6;k++)
				alpha[j] -= M[k*3+j]*x[k+7]; // 由极大值原理确定推力方向
		}
		norma = enorm(3, alpha);
		for(int j=0;j<3;j++)
			alpha[j] /= norma;
		rou = 1.0-(Ispg0NU*norma/x[6] + x[13])/Out[1]; // 开关函数
		// u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // 对数型同伦指标下的推力
		if (rou>0)
			u = 0;
		else
			u = 1;
		fprintf(ufid, "%.6f, %.6f\n", TimeNode[i], u);
	}
	fclose(ufid);
}