#include"OrbitFun.h"
//namespace OrbitFun{
//#define DPI 3.1415926535897932384626433832795
//#define D2PI 2.0*PI

// *             �廪��ѧ���캽��ѧԺ����ѧ������о��Ҳ�ʿ��               *
// *                ������(jiangfh04@mails.thu.edu.cn)                      *
// *                 ����޸�: 2009.2.24                                    *


// ���ݾ�������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
//������� coe[5],mu:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14

//�������RV:
//   λ����ʸ��R����X��Y��Z���ٶ���ʸ��V�ķ���VX��VY��VZ��6ά����
//   ��λ���ס���ÿ��.
//�ɹ�flag����1������0
void coe2rv(int& flag, double* rv, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<=0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)
		return;
	if((coe[1]*cos(coe[5]))<-1.0)
	{
//		cout<<"�����ܴﵽ��˫�����."<<endl;		
		return;
	}

	double p=coe[0]*fabs(1.0-coe[1]*coe[1]);//��ͨ��
	if(coe[1]==1.0)//����������߹��,����Դ�.
		p=2.0*coe[0];	

	double sini, cosi, sinO, cosO, sino, coso;
	sini=sin(coe[2]);
	cosi=cos(coe[2]);
	sinO=sin(coe[3]);
	cosO=cos(coe[3]);
	sino=sin(coe[4]);
	coso=cos(coe[4]);

	//���ƽ�淨��λʸ��,���Ƕ�����λʸ��
	double HVector[3]={sini*sinO, -sini*cosO, cosi};
      
	//ƫ���ʵ�λʸ��,���Laplaceʸ��
	double PVector[3]={cosO*coso-sinO*sino*cosi, sinO*coso+cosO*sino*cosi, sino*sini};

	//��ͨ������λʸ��,PVector,QVector,HVector������������ϵ
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r=0.0;
	if((coe[1]*cos(coe[5]))+1.0<=0.0)
	{
//		cout<<"�����˫�������������Զ��."<<endl;
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

// E2f ����ƫ����Ǻ�ƫ������������
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH,��r=a(ecoshH-1)
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������f:
// 	�����ǣ����ȣ�
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
double E2f(int& flag, double E, double e)
{
	if(e<0.0) {flag=0;return E;}

	double f=0.0;
    if(e>=0.0&&e<1.0)//Բ����Բ���
	{
		double E0=fmod(E, D2PI);
		if(E0>DPI)
			E0-=D2PI;
		if(E0<-DPI)
			E0+=D2PI;
		f=2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*E0));
		f=f+E-E0;
	}
	else if(e>1.0)//˫���߹��
		f=2.0*atan(sqrt((e+1.0)/(e-1.0))*tanh(0.5*E));
	else// (abs(e-1.0)<epsilon)�����߹��
	{
		f=E;
//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽������ǵ�ֵ��ΪE."<<endl;
	}
	flag=1;
	return f;	
}


// E2M ����ƫ����Ǻ�ƫ������ƽ�����
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������M:
// 	ƽ�����(����),����˫�����,ָN
double E2M(int& flag, double E, double e)
{
	if(e<0.0){flag=0;return E;}
	double M=0.0;
	if(e>=0.0&&e<1.0)//Բ����Բ���
	{
		double E0=fmod(E, D2PI);    
		M=E0-e*sin(E0);
		M=M+E-E0;
	}
	else if(e>1.0)//˫���߹��
		M=e*sinh(E)-E;   
	else//(abs(e-1.0)<epsilon)�����߹��
	{
		M=E;
//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽�ƽ�����ֵ��ΪE."<<endl;
	}
	flag=1;
	return M;
}

// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
//�������f0,t,a,e,mu,MaxIter,epsilon:
//   f0:��ʼ������,��λ������
//   t:����ʱ��,��λ����
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   MaxIter:����������
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������ft:
// 	����������(����)
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if(mu<=0.0||MaxIter<1||a<=0.0||e<0.0){flag=0;return f0;}

	double ft=0.0;
	if((e>=0.0&&e<1.0)||(e>1.0))//Բ,��Բ,˫�����
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
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		if((f0<-DPI)||(f0>DPI))
		{
//			cout<<"���������߹������ʼ������Ӧ��-180��180��֮��."<<endl;
			flag=0;return f0;
		}
		else if(f0>DPI||f0<-DPI)
			ft=f0;
		else
		{
			double B=0.75*sqrt(2.0*mu/(a*a*a))*dt+0.5*tan(0.5*f0)*((tan(0.5*f0))*(tan(0.5*f0))+3.0);
			double B1B=B+sqrt(1.0+B*B);
			double tanv=0.0;
			if(fabs(dt)<D2PI*sqrt((a*a*a)/mu)/1000.0)//�ƽ�ʱ��ΪС�������
			{
				double A=pow(B1B, 2.0/3.0);
				tanv=2.0*A*B/(1.0+(1.0+A)*A);
			}
			else//����С�������
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


// f0ft2dt ���ݳ�ʼ�����Ǻ��������������ݻ�ʱ��
//�������f0,ft,a,e,mu,AbsTol:
//   f0:��ʼ������,��λ������
//   ft:�ն�������,��λ������
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������t:
//   t:����ƽ�ʱ��,��λ����
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
//			cout<<"�����ܴﵽ��˫���������߹��"<<endl;			
			return 0.0;
		}
		else if(f0<-maxangle||f0>maxangle||ft<-maxangle||ft>maxangle)
		{
//			cout<<"����ʱ����Ѽ���׼ȷ,����Ϊ����,�ڴ���Ϊ1.0e100.";
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


// f2E ���������Ǻ�ƫ������ƫ�����
//�������f, e:
//   f:������,��λ������
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������E:
// 	ƫ����ǣ����ȣ�,����˫�����,ָH,��r=a(ecoshH-1)
double f2E(int& flag, double f, double e)
{
	if(e<0.0){flag=0;return f;}

	double E=0.0;
	if(e>=0.0&&e<1.0)//Բ����Բ���
	{
		double f0=fmod(f, D2PI);
		if(f0>DPI)
			f0-=D2PI;
		if(f0<-DPI)
			f0+=D2PI;
		E=2.0*atan(sqrt((1.0-e)/(1.0+e))*tan(0.5*f0));
		E+=f-f0;
	}
	else if(e>1.0)//˫���߹��
	{
		if(f>DPI-acos(1.0/e)||f<-DPI+acos(1.0/e))
		{
//			cout<<"�����ܴﵽ��˫�����."<<endl;
			flag=0;
			return f;
		}
		else
			E=2.0*atanh(sqrt((e-1.0)/(1.0+e))*tan(0.5*f));
	}
	else//if(abs(e-1.0)<epsilon)�����߹��
	{
		E=f;
//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��Ϊf."<<endl;
	}
	flag=1;
	return E;
}


// Inertia2LVLH �������ǵľ���λ��R0���ٶ�V0�����ǵľ���λ��R1���ٶ�V1����������ǹ������ϵ�е����λ�ú��ٶ�
//ע:R0,V0,R1,V1�������ϵ��ͶӰ,r,v�����ǹ������ϵͶӰ
//�������R0[2], V0[2], R1[2], V1[2], epsilon:
//   R0[2], R1[2]:���Ǻʹ�����Ե��ĵ�λ����ʸ������X��Y��Z,��λ����
//   V0[2], V1[2]:���Ǻʹ�����Ե��ĵ��ٶ���ʸ������VX��VY��VZ,��λ����ÿ��.
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������rv[5]:���������ǹ������ϵ�е����λ�ú�����ٶ���ʸ������,��λ����,��ÿ��
void Inertia2LVLH(int& flag, double* rv, const double* RV0, const double* RV1)
{
	int i;
	
	double R0[3]={RV0[0], RV0[1], RV0[2]};
	double V0[3]={RV0[3], RV0[4], RV0[5]};
	double R1[3]={RV1[0], RV1[1], RV1[2]};
	double V1[3]={RV1[3], RV1[4], RV1[5]};
	double radius=V_Norm2(R0,3);//����
	if(radius<=0.0){flag=0;return;}
	double UnitR[3];
	for(i=0;i<3;i++) UnitR[i]=R0[i]/radius;//����λʸ��    
	double VectorH[3];
	V_Cross(VectorH, R0, V0);
	double h=V_Norm2(VectorH,3);//�Ƕ���ֵ
	if(h<=0.0){flag=0;return;}
	double UnitH[3];
	for(i=0;i<3;i++) UnitH[i]=VectorH[i]/h;//����淨��λʸ��
	
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

// LVLH2Inertia �������ǵľ���λ��R0���ٶ�V0������������ǹ������ϵ�е�λ��r���ٶ�v����ǵľ���λ�ú��ٶ�
//ע:RV0,RV1�����ϵͶӰ;rv�����ǹ������ϵͶӰ
//�������RV0[5],rv[5],epsilon:
//   RV0[5]:������Ե��ĵ�λ�ú��ٶ���ʸ������X��Y��Z,��λ����,VX��VY��VZ,��λ����ÿ��
//   rv[5]:����������ǹ������ϵ�е�λ�ú��ٶ���ʸ������X��Y��Z,��λ����,VX��VY��VZ,��λ����ÿ��
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������RV1[5]:������Ե��ĵ�λ�ú��ٶ���ʸ������X��Y��Z,��λ����,VX��VY��VZ,��λ����ÿ��
void LVLH2Inertia(int& flag, double* RV1, const double* RV0, const double* rv, double epsilon)
{
	int i;
	
	double R0[3]={RV0[0], RV0[1], RV0[2]};
	double V0[3]={RV0[3], RV0[4], RV0[5]};
	double r[3]={rv[0], rv[1], rv[2]};
	double v[3]={rv[3], rv[4], rv[5]};
	double radius=V_Norm2(R0, 3);//����
	if(radius<=0.0){flag=0;return;}
	double UnitR[3];
	for(i=0;i<3;i++) UnitR[i]=R0[i]/radius;//����λʸ��    
	double VectorH[3];
	V_Cross(VectorH, R0, V0);
	double h=V_Norm2(VectorH,3);//�Ƕ���ֵ
	if(h<=0.0){flag=0;return;}
	double UnitH[3];
	for(i=0;i<3;i++) UnitH[i]=VectorH[i]/h;//����淨��λʸ��
	
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


// M2E ����ƽ����Ǻ�ƫ������ƫ�����
//�������M, e, MaxIter, epsilon:
//   M:ƽ�����,��λ������,����˫�����,ָN
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   MaxIter:����������,Ĭ��Ϊ60
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������E:
// 	ƫ����ǣ����ȣ�����˫�������ָH
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if(epsilon<=0.0||MaxIter<1||e<0.0){flag=0;return M;}

	//������������Solar System Dynamics��Chapter2,Carl D.Murray and Stanley F.Dermott��
	double E=0.0, Minus=0.0, DeMinus=0.0, DeDeMinus=0.0, DeDeDeMinus=0.0, Delta1=0.0, Delta2=0.0, Delta3=0.0;
	int N=0;
	if(e>=0.0&&e<1.0)//Բ����Բ���
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
	else if(e>1.0)//˫���߹��
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
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		E=M;
//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��ΪM."<<endl;
	}
	if(((e>=0.0&&e<1.0)||(e>1.0))&&fabs(Delta3)>=5.0*epsilon&&N>=MaxIter)
	{
//		cout<<"����������,�뽵�;���epsilon�����ӵ�����������."<<endl;
		flag=0;
		return M;
	}
	flag=1;
	return E;
}



// rv2coe ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ����󾭵�������
//�������RV[5], mu, epsilon:
//   RV:λ����ʸ������X��Y��Z,��λ����,���ٶ���ʸ������VX��VY��VZ,��λ����ÿ��.
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������coe[6]:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�q,p=2q.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,����ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵȡΪomega+f,��γ�ȷ���
void rv2coe(int& flag, double* coe, const double* RV, double mu)
{
	int i;
	flag=0;
	if(mu<=0.0)
		return;

	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R,3);//����
	double velocity=V_Norm2(V,3);//�ٶ�
	if(radius<=0.0||velocity<=0.0)
		return;
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//����λʸ��    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector,unitR,unitV);
	double h=radius*velocity*V_Norm2(hvector,3);//�Ƕ���ֵ
	if(h<=0.0)
		return;
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector,3);//����淨��λʸ��
	//ƫ����ʸ��
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];
	coe[1]=V_Norm2(evector,3);//ƫ����
	double p=h*h/mu;
	if(coe[1]==1.0)
		coe[0]=0.5*p;//�����߹���Ľ��Ǿ�
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//�볤��
	bool judge=(coe[1]>0.0);
	double unite[3]={0.0};
	if(judge)
		for(i=0;i<3;i++) unite[i]=evector[i]/coe[1];//ƫ���ʵ�λʸ��
	coe[2]=acos(unith[2]);//������

	double unitN[3]={-unith[1], unith[0], 0.0};//����ʸ��,δ��һ��

	double temp[3];

	if(V_Norm2(unitN,3)==0.0)
	{
		coe[3]=0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if(!judge)
		{
			coe[4]=0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5]=atan2(unitR[1]*unith[2],unitR[0]);//������
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4]=atan2(unite[1]*unith[2],unite[0]); //���ǵ����       
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
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4]=atan2(V_Dot(unith,temp,3), V_Dot(unite,unitN,3));
			coe[5]=coe[5]-coe[4];
		}
	}
	//ת����[0,2pi)��
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
	am=(V_Norm2(rv0,3)+V_Norm2(rv1,3)+V_Norm2(vt0,3))/4.0;//��С��Բת�ƹ���볤��
	Nmax=(int)floor(t/(D2PI*sqrt(am*am*am/GM)));//�����ܵ�ת��Ȧ��
	for(n=0;n<=Nmax;n++)//Ȧ����0������������٣����������ٶ���С�Ľ��
	{
		for(int j=0;j<2;j++)//���ڶ�Ȧ���⣬�����ҷ�֦������
		{
			lambert(v0, v1, a, e, rv0, rv1, t, unith, flag0, GM, 0, n, j);
			if(flag0!=1||e>=1.0) continue;//�ų�˫����������
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

// coe2ee ���ݾ������������������������Թ�����180������
//������� coe[5],mu:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������ee[6]:
// 	ee[0]:��ͨ��p���ף�.
// 	ee[1]:e*cos(omega+OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	ee[2]:e*sin(omega+OMEGA),��g��ʾ
//	ee[3]:tan(i/2)*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	ee(5:tan(i/2)*sin(OMEGA),��k��ʾ
//	ee[5]:omega+OMEGA+f,��L��ʾ
void coe2ee(int&flag, double* ee, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)	
		return;	
	if((coe[1]*cos(coe[5]))<-1.0)
	//		cout<<"�����ܴﵽ��˫�����."<<endl;
		return;
	
	ee[0]=coe[0]*fabs(1.0-coe[1]*coe[1]);//��ͨ��p
	if(coe[1]==1.0)//����������߹��,����Դ�.
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


// coe2mee ���ݾ�������������������������
//������� coe[5],mu:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������[mee[6],orbtype]:
// 	mee[0]:��ͨ��p���ף�.
// 	mee[1]:e*cos(omega+I*OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	mee[2]:e*sin(omega+I*OMEGA),��g��ʾ
//	mee[3]:(tan(i/2))^I*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	mee(5:(tan(i/2))^I*sin(OMEGA),��k��ʾ
//	mee[5]:omega+I*OMEGA+f,��L��ʾ
//   orbtype:˳�й��1,���й��-1
void coe2mee(int&flag, double* mee, int& orbtype, const double* coe, double mu)
{
	flag=0;
	if(mu<=0.0||coe[0]<0.0||coe[1]<0.0||coe[2]<0.0||coe[2]>DPI)	
		return;	
	if((coe[1]*cos(coe[5]))<-1.0)	
//		cout<<"�����ܴﵽ��˫�����."<<endl;
		return;

	orbtype=1;
	if(coe[2]>DPI/2.0)
		orbtype=-1;

	mee[0]=coe[0]*fabs(1.0-coe[1]*coe[1]);//��ͨ��p
	if(coe[1]==1.0)//����������߹��,����Դ�.
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

// dee_dt ������������ee���ؾ��򡢺��򡢷�����ٶ�Ad��������ʱ��ı仯��
//������� Ad, ee[5],mu:
//   Ad[0]:������ٶȴ�С��
//   Ad[1]:������ٶȴ�С��
//   Ad[0]:������ٶȴ�С����λm/s^2
// 	ee[0]:��ͨ��p���ף�.
// 	ee[1]:e*cos(omega+OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	ee[2]:e*sin(omega+OMEGA),��g��ʾ
//	ee[3]:tan(i/2)*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	ee(5:tan(i/2)*sin(OMEGA),��k��ʾ
//	ee[5]:omega+OMEGA+f,��L��ʾ
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//�������dee[6]:�ֱ���ee[6]ÿ��������Ӧ����ʱ��仯��
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


// dmee_dt ����������������mee���ؾ��򡢺��򡢷�����ٶ�Ad��������ʱ��ı仯��
//������� Ad, mee[5],orbtype,mu:
//   Ad[0]:������ٶȴ�С��
//   Ad[1]:������ٶȴ�С��
//   Ad[0]:������ٶȴ�С����λm/s^2
// 	mee[0]:��ͨ��p���ף�.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	mee[2]:e*sin(omega+orbtype*OMEGA),��g��ʾ
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	mee[4]:(tan(i/2))^orbtype*sin(OMEGA),��k��ʾ
//	mee[5]:omega+orbtype*OMEGA+f,��L��ʾ
//   orbtype:�����ʶ,˳�й��1,���й��-1,Ĭ��Ϊ˳�й��
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//�������dmee[6]:�ֱ���mee[6]ÿ��������Ӧ����ʱ��仯��
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

// ee2coe ������������󾭵�������
//������� ee[5],mu,epsilon:
// 	ee[0]:��ͨ��p���ף�.
// 	ee[1]:e*cos(omega+OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	ee[2]:e*sin(omega+OMEGA),��g��ʾ
//	ee[3]:tan(i/2)*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	ee(5:tan(i/2)*sin(OMEGA),��k��ʾ
//	ee[5]:omega+OMEGA+f,��L��ʾ
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������coe[6]:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�q,p=2q.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,����ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵȡΪomega+f,��γ�ȷ���
void ee2coe(int& flag, double* coe, const double* ee, double mu)
{
	flag=0;
	if(mu<=0.0||ee[0]<=0.0)
		return;

	double p=ee[0], f=ee[1], g=ee[2], h=ee[3], k=ee[4], L=ee[5];
	coe[1]=sqrt(f*f+g*g);
	if(coe[1]==1.0)
		coe[0]=0.5*p;//�����߹���Ľ��Ǿ�
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//�볤��
	double temp=sqrt(h*h+k*k);
	coe[2]=2.0*atan(temp);
	if(temp<=0.0)
	{
		coe[3]=0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5]=L;//������
		}
		else
		{
			coe[4]=atan2(g,f); //���ǵ����       
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
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			coe[4]=atan2(g*h-f*k,f*h+g*k);
			coe[5]=coe[5]-coe[4];
		}
	}
	//ת����[0,2pi)��
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


// ee2rv ������������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
//������� ee[5],mu:
// 	ee[0]:��ͨ��p���ף�.
// 	ee[1]:e*cos(omega+OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	ee[2]:e*sin(omega+OMEGA),��g��ʾ
//	ee[3]:tan(i/2)*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	ee(5:tan(i/2)*sin(OMEGA),��k��ʾ
//	ee[5]:omega+OMEGA+f,��L��ʾ
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//�������RV:
//   λ����ʸ��R����X��Y��Z���ٶ���ʸ��V�ķ���VX��VY��VZ��ɵ�һ��6ά������,
//   ��λ���ס���ÿ��.
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

// mee2coe �������������������󾭵�������
//������� mee[5],orbtype,mu,epsilon:
// 	mee[0]:��ͨ��p���ף�.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	mee[2]:e*sin(omega+orbtype*OMEGA),��g��ʾ
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),��k��ʾ
//	mee[5]:omega+orbtype*OMEGA+f,��L��ʾ
//   orbtype:�����ʶ,˳�й��1,���й��-1,Ĭ��Ϊ˳�й��
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������coe[6]:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�q,p=2q.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,����ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,����ֵȡΪomega+f,��γ�ȷ���
void mee2coe(int&flag, double* coe, const double* mee, int orbtype, double mu)
{
	flag=0;
	if(mu<=0.0||abs(orbtype)!=1||mee[0]<=0.0)
		return;

	double p=mee[0], f=mee[1], g=mee[2], h=mee[3], k=mee[4], L=mee[5];
	coe[1]=sqrt(f*f+g*g);
	if(coe[1]==1.0)
		coe[0]=0.5*p;//�����߹���Ľ��Ǿ�
	else
		coe[0]=p/(fabs(1.0-coe[1]*coe[1]));//�볤��
	double temp=sqrt(h*h+k*k);
	if(temp>1.0)
	{
//		cout<<"��4��5���������ƽ���Ͳ��ܴ���1."<<endl;
		return;
	}
	if(orbtype==1)
		coe[2]=2.0*atan(temp);
	else
		coe[2]=2.0*(DPI/2.0-atan(temp));
	if(temp<=0.0)
	{
		coe[3]=0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if(coe[1]<=0.0)
		{
			coe[4]=0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5]=L;//������
		}
		else
		{
			coe[4]=atan2(g,f); //���ǵ����       
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
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			coe[4]=atan2(g*h-orbtype*f*k,f*h+orbtype*g*k);
			coe[5]=coe[5]-coe[4];
		}
	}
	//ת����[0,2pi)��
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


// mee2rv ����������������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
//������� mee[5],orbtype,mu:
// 	mee[0]:��ͨ��p���ף�.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	mee[2]:e*sin(omega+orbtype*OMEGA),��g��ʾ
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),��k��ʾ
//	mee[5]:omega+orbtype*OMEGA+f,��L��ʾ
//   orbtype:�����ʶ,˳�й��1,���й��-1,Ĭ��Ϊ˳�й��
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//�������RV:
//   λ����ʸ��R����X��Y��Z���ٶ���ʸ��V�ķ���VX��VY��VZ,
//   ��λ���ס���ÿ��.
void mee2rv(int&flag, double* rv, const double* mee, int orbtype, double mu)
{
	flag=0;
	if(mu<=0.0||abs(orbtype)!=1||mee[0]<=0.0)
		return;	

	double p=mee[0], f=mee[1], g=mee[2], h=mee[3], k=mee[4], L=mee[5];
	double temp=sqrt(h*h+k*k);
	if(temp>1.0)
	{
//		cout<<"��4��5���������ƽ���Ͳ��ܴ���1."<<endl;
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

// rv2ee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ������������������Թ�����180��ʱ����
//�������RV[5], mu, epsilon:
//   RV:λ����ʸ������X��Y��Z,��λ����,���ٶ���ʸ������VX��VY��VZ,��λ����ÿ��.
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������ee[6]:
// 	ee[0]:��ͨ��p���ף�.
// 	ee[1]:e*cos(omega+OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	ee[2]:e*sin(omega+OMEGA),��g��ʾ
//	ee[3]:tan(i/2)*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	ee(5:tan(i/2)*sin(OMEGA),��k��ʾ
//	ee[5]:omega+OMEGA+f,��L��ʾ
void rv2ee(int&flag, double* ee, const double* RV, double mu)
{
	flag=0;
	if(mu<=0.0)
		return;
	int i;
	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R, 3);//����
	double velocity=V_Norm2(V, 3);//�ٶ�	
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//����λʸ��    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h=radius*velocity*V_Norm2(hvector, 3);//�Ƕ���ֵ
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector, 3);//����淨��λʸ��
	
	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//ƫ����ʸ��	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];	
	//���ܾ�����ʸ��,ģΪe
	double qvector[3];
	V_Cross(qvector,unith,evector);
	//�����һ��ʸ��
	double unitA[3];
	for(i=0;i<3;i++) unitA[i]=h/mu*V[i]-qvector[i];
	//������ϵʽ
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
	//�Ƶ��ó�
	//evector[0]+qvector[1]=(1+cos(i))*e*cos(omega+OMEGA);
	//evector[1]-qvector[0]=(1+cos(i))*e*sin(omega+OMEGA);
	//unith[0]=(1+cos(i))*tan(i/2)*sin(OMEGA);
	//unith[1]=-(1+cos(i))*tan(i/2)*cos(OMEGA);
	// unitR[0]+unitA[1]=(1+cos(i))*cos(omega+OMEGA+f);
	// unitR[1]-unitA[0]=(1+cos(i))*sin(omega+OMEGA+f);

	ee[0]=h*h/mu;//p

	if(unith[2]+1.0<=0.0)
	{
//		cout<<"�����ǽӽ�180�ȣ����ʺ�����������������. "<<endl;
		return;
	}
	double cosiadd1=1.0+unith[2];
	ee[1]=(evector[0]+qvector[1])/cosiadd1;//f
	ee[2]=(evector[1]-qvector[0])/cosiadd1;//g
	ee[3]=-unith[1]/cosiadd1;//h
	ee[4]=unith[0]/cosiadd1;//k
	ee[5]=atan2(unitR[1]-unitA[0],unitR[0]+unitA[1]);//L
	//ת����[0,2pi)��
	ee[5]=fmod(ee[5], D2PI);
	if(ee[5]<0.0)
		ee[5]+=D2PI;
	flag=1;
	return;
}


// rv2mee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ�������������������
//�������RV[5], mu, epsilon:
//   RV:λ����ʸ������X��Y��Z,��λ����,���ٶ���ʸ������VX��VY��VZ,��λ����ÿ��.
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������[mee[6],orbtype]:
// 	mee[0]:��ͨ��p���ף�.
// 	mee[1]:e*cos(omega+orbtype*OMEGA),��f(��Ҫ�������ǻ���)��ʾ
// 	mee[2]:e*sin(omega+orbtype*OMEGA),��g��ʾ
//	mee[3]:(tan(i/2))^orbtype*cos(OMEGA),��h(��Ҫ��Ƕ�����С����)��ʾ
//	mee(5:(tan(i/2))^orbtype*sin(OMEGA),��k��ʾ
//	mee[5]:omega+orbtype*OMEGA+f,��L��ʾ
//   orbtype:˳�й��1,���й��-1
void rv2mee(int&flag, double* mee, int& orbtype, const double* RV, double mu, double epsilon)
{
	flag=0;
	if(mu<=0.0)
		return;
	int i;

	double R[3]={RV[0], RV[1], RV[2]};
	double V[3]={RV[3], RV[4], RV[5]};
	double radius=V_Norm2(R, 3);//����
	double velocity=V_Norm2(V, 3);//�ٶ�
	assert(radius>=epsilon&&velocity>=epsilon);
	double unitR[3];
	for(i=0;i<3;i++) unitR[i]=R[i]/radius;//����λʸ��    
	double unitV[3];
	for(i=0;i<3;i++) unitV[i]=V[i]/velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector,unitR,unitV);
	double h=radius*velocity*V_Norm2(hvector, 3);//�Ƕ���ֵ
	assert(h>=epsilon);
	double unith[3];
	for(i=0;i<3;i++) unith[i]=hvector[i]/V_Norm2(hvector, 3);//����淨��λʸ��
	
	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//ƫ����ʸ��	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for(i=0;i<3;i++) evector[i]=(velocity*h/mu)*evector[i]-unitR[i];	
	//���ܾ�����ʸ��,ģΪe
	double qvector[3];
	V_Cross(qvector,unith,evector);
	//�����һ��ʸ��
	double unitA[3];
	for(i=0;i<3;i++) unitA[i]=h/mu*V[i]-qvector[i];
	//������ϵʽ
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
	//�Ƶ��ó�
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
	//ת����[0,2pi)��
	mee[5]=fmod(mee[5], D2PI);
	if(mee[5]<0.0)
		mee[5]+=D2PI;
	flag=1;
	return;
}

//}//end namespace OrbitFun
