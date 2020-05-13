#include "ExternHeader.h"
#include "Constant.h"

int Derivative(double t, const double* x, double* dx, const double* dfpara)
{
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // ͬ�ײ���
	double lam0 = dfpara[1];

	// ������������λʸ��
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	// ����������ֵ��С
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // ������ͬ��ָ���µ�����

	double r = V_Norm2(rv, 3);
	for (int i=0;i<3;++i)
	{
		dx[i] = x[i+3]; // 3��λ�÷����ĵ�����3���ٶȷ���
		dx[i+3] = -muNU/(r*r*r)*rv[i] + u*TmaxNU/m*alpha[i];
		
		dx[i+7] = muNU/(r*r*r)*costate[i+3] - 3*muNU*V_Dot(rv, &costate[3], 3)/pow(r, 5)*rv[i];
		dx[i+10] = -costate[i];
	}
	dx[6] = -u*TmaxNU/Ispg0NU;
	dx[13] = -u*TmaxNU/(m*m)*norm_lambdaV;
	return 1;
}


//����Derivative()�õ�14�����ֱ�����ʱ��仯�ĵ�����ͨ��ode45()���и���Э̬��ֵ����µĻ��֣��õ����Ƶ��ն�״̬����
//����ʱ���״̬�̶���ȼ�����������з���
//sfpara
//[0-6],��ʼλ���ٶȺ�����
//[7-12],ĩ��λ���ٶ�
//[13-14],ת��ʱ���epsilon
//[15],�����Ϣ��ʶ��0��ʾ�������1��ʾ���
//x:[0],lam0;[1-7]��ʼЭ̬
//nΪ����ά����fvecΪ�вiflagΪ0ʱ��ʾ�Ϸ����ú�����֪������в�
//����ɹ�ʱ������0�����0��ֵ�����ɹ�ʱ����С��0��ֵ
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

	// FILE *fid = NULL;//fopen("temp.txt","w");//����趨��Ч�ļ�·���������Ҫ�ر��ļ�
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
	double sfpara[16] = {0.0}; // 0~6-7����ʼ״ֵ̬��7~13-7���ն�״ֵ̬��14-ͬ�ײ�����15-ShootFunFuelOpt()�������Ϣ��־
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
				x[j]=(double)rand()/RAND_MAX-0.5; // ���������б�����ֵ
		x[7]+=0.5; // ����Э̬����Ϊ������ĩ��Ϊ0����˳�ֵ��Ϊ��
		x[0]+=0.5; // lambda0����һ�����ӣ���Ϊȡ��
		info = hybrd1(ComputeFvec, n, x, fvec, sfpara, wa, xtol, 20, 500); // info-hybrd1()�������־
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			int _j = ComputeFvec(n, x, fvec, 1, sfpara); // ����ͬ�׼���õ���Э̬��ֵ���������һ�εĻ�����⣨ֱ���������С��ͬ�ײ����������ư����ƵĽ����
			if(fvec[0]>Out[0]) // ʣ������Ϊ�����������о���Ҫ��ֹͣ
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