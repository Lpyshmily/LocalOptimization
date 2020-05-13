#include "OrbitFun.h"


// E2f ����ƫ����Ǻ�ƫ������������
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH,��r=a(ecoshH-1)
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������f:
// 	�����ǣ����ȣ�
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
double E2f(int& flag, double E, double e)
{
	if(e<0.0) {flag=0; return E;}

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
//�������M:
// 	ƽ�����(����),����˫�����,ָN
double E2M(int& flag, double E, double e)
{
	if(e<0.0){flag=0; return E;}
	double M=0.0;
	if(e>=0.0&&e<1.0)//Բ����Բ���
	{
		double E0=fmod(E, D2PI);    
		M=E0-e*sin(E0);
		M=M+E-E0;
	}
	else if(e>1.0)//˫���߹��
		M=e*sinh(E)-E;   
	else // (abs(e-1.0)<epsilon)�����߹��
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
	if(mu<=0.0||MaxIter<1||a<=0.0||e<0.0) {flag=0; return f0;}

	double ft=0.0;
	if((e>=0.0&&e<1.0)||(e>1.0))//Բ,��Բ,˫�����
	{
		double E = f2E(flag,f0,e);
		if(flag==0) return f0;
		double M = E2M(flag,E,e);
		if(flag==0) return f0;
		M += sqrt(mu/(a*a*a))*dt;
		E = M2E(flag,M,e,MaxIter,epsilon);
		if(flag==0) return f0;
		ft = E2f(flag,E,e);
		if(flag==0) return f0;
	}
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		if((f0<-DPI)||(f0>DPI))
		{
//			cout<<"���������߹������ʼ������Ӧ��-180��180��֮��."<<endl;
			flag=0; return f0;
		}
		else if(f0>DPI||f0<-DPI)
			ft = f0;
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
//�������flag,f0,ft,a,e,mu,AbsTol:
//   f0:��ʼ������,��λ������
//   ft:�ն�������,��λ������
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
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
		double maxangle=DPI-acos(1.0/e); // ˫����/�����������������ǣ������ߣ�
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
	else // if(fabs(e-1.0)<epsilon)
	{
		double B1=tan(0.5*f0)*((tan(0.5*f0))*(tan(0.5*f0))+3.0); // 2B
		double B2=tan(0.5*ft)*((tan(0.5*ft))*(tan(0.5*ft))+3.0); // 2B
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
	if(e<0.0){flag=0; return f;}

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
		while(fabs(Delta3)>=epsilon && N<MaxIter)
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
	if(((e>=0.0&&e<1.0)||(e>1.0)) && fabs(Delta3)>=5.0*epsilon && N>=MaxIter)
	{
//		cout<<"����������,�뽵�;���epsilon�����ӵ�����������."<<endl;
		flag=0;
		return M;
	}
	flag=1;
	return E;
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