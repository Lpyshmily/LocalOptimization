#include "ExternHeader.h"
#include "GA_FOP.h"

// ����ͬ��
// ����14�����ֱ���x���������Ӧ�ĵ���dx
// dfpara:[0],epsi;[1],lam0
// ��һ������tֻ��Ϊ��ode45���ã��ڳ����в�������
int GA_derivative(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // ͬ�ײ���
	double lam0 = dfpara[1];

	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // ������ͬ��ָ���µ�����

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
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

// bang-bang����
// ����14�����ֱ���x���������Ӧ�ĵ���dx
// dfpara:[0],lam0
int GA_derivative_type1(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // ͬ�ײ���
	double lam0 = dfpara[1];

	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	if (rou > 0)
		u = 0.0;
	else
		u = 1.0;
	

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
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

// ����ͬ��
int GA_derivative_type2(double t, const double* x, double* dx, const double* dfpara)
{
	int i;
	
	double rv[6] = {0.0}, costate[7] = {0.0}, alpha[3] = {0.0};
	double m;
	V_Copy(rv, x, 6);
	m = x[6];
	V_Copy(costate, &x[7], 7);
	double epsi = dfpara[0]; // ͬ�ײ���
	double lam0 = dfpara[1];

	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/2.0/epsi;

	double r = V_Norm2(rv, 3);
	for (i=0;i<3;++i)
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

// ����ͬ��
// ������ܶٺ�����ֵ
double GA_Hamilton(const double* x, double lam0, double epsi)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	double u = 2.0*epsi/(rou+2.0*epsi+sqrt(rou*rou+4.0*epsi*epsi)); // ������ͬ��ָ���µ�����
	// ����״̬�����ĵ���
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

// bang-bang����
// ������ܶٺ�����ֵ
double GA_Hamilton_type1(const double* x, double lam0)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	if (rou > 0)
		u = 0.0;
	else
		u = 1.0;
	// ����״̬�����ĵ���
	double dfpara = lam0;
	double dx[14] = {0.0};
	int flag = GA_derivative_type1(0, x, dx, &dfpara);

	double H = 0.0;
	H = V_Dot(costate, dx, 7);
	H = H + lam0*TmaxNU/Ispg0NU*u;
	return H;
}

// ����ͬ��
double GA_Hamilton_type2(const double* x, double lam0, double epsi)
{
	double m = x[6];
	double costate[7] = {0.0}, alpha[3] = {0.0};
	V_Copy(costate, &x[7], 7);
	
	// ������������λʸ����������ֵ��С
	double norm_lambdaV = V_Norm2(&costate[3], 3); // �ٶ�Э̬�ķ���
	for (int i=0;i<3;++i)
		alpha[i] = -costate[i+3]/norm_lambdaV;
	double rou, u;
	rou = 1.0 - (Ispg0NU*norm_lambdaV/m + costate[6])/lam0; // ���غ���
	if (rou > epsi)
		u = 0.0;
	else if (rou < -epsi)
		u = 1.0;
	else
		u = 0.5 - rou/2.0/epsi;
	// ����״̬�����ĵ���
	double dfpara[2];
	dfpara[0] = epsi;
	dfpara[1] = lam0;
	double dx[14] = {0.0};
	int flag = GA_derivative_type2(0, x, dx, dfpara);

	double H = 0.0;
	H = V_Dot(costate, dx, 7);
	H = H + (lam0*TmaxNU/Ispg0NU*( u-epsi*log(u*(1.0-u)) ));
	return H;
}

//����GA_derivative()�õ�14�����ֱ�����ʱ��仯�ĵ�����ͨ��ode45()���и���Э̬��ֵ����µĻ��֣��õ����Ƶ��ն�״̬����
//����ʱ���״̬�̶���ȼ�����������з���
//sfpara
//[0-6],��ʼλ���ٶȺ�����
//[7-12],ĩ��λ���ٶ�
//[13-14],ת��ʱ���epsilon
//[15],�����Ϣ��ʶ��0��ʾ�������1��ʾ���
//[16-21],������MJD_MARSʱ�̵�λ���ٶ�
//x:[0],lam0;[1-7]��ʼЭ̬;[8-11]��ʽԼ������;[12]����ʽԼ������;[13-15]���������ٶ�����dvg;[16]��������ʱ��
//nΪ����ά����fvecΪ�вiflagΪ0ʱ��ʾ�Ϸ����ú�����֪������в�
//����ɹ�ʱ������0�����0��ֵ�����ɹ�ʱ����С��0��ֵ
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
	// ��ȡMJD_MARSʱ�̻���λ��
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

	// FILE *fid = NULL;//fopen("temp.txt","w");//����趨��Ч�ļ�·���������Ҫ�ر��ļ�
	// FILE *fid1 = fopen("temp1.txt", "w");
	// FILE *fid2 = fopen("temp2.txt", "w");
	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	int flag,NumPoint;
	double work[140]={0.0};
	flag = ode45(GA_derivative, x0, dfpara,  0.0, tm*TOF*86400/TUnit, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid1);
	// ����һ����ƫ��
	// ��������λ��Լ��,7-9
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
	// �����ɳ�����,11
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

	// tm-ʱ���ٶ�Э̬,12-14
	for (int i=0;i<3;++i)
		fvec[12+i] = x0[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	double H1 =  GA_Hamilton(x0, lam0, epsi);
	// �����µ�״̬������Э̬����
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // �ٶ�
		x0[7+i] = x0[7+i] - chi[i]; // λ��Э̬
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // �ٶ�Э̬
	}
	double H2 = GA_Hamilton(x0, lam0, epsi);
	// ��̬����ƫ��,15
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - kappa*c;


	// �ڶ��׶εĻ���
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

//sfpara
//[0-6],��ʼλ���ٶȺ�����
//[7-12],ĩ��λ���ٶ�
//[13],ת��ʱ��
//[14],�����Ϣ��ʶ��0��ʾ�������1��ʾ���
//[15-20],������MJD_MARSʱ�̵�λ���ٶ�
int GA_fvec_type1(int n, const double *x, double *fvec, int iflag, const double* sfpara)
{
	if(iflag==0)
	{
		return 0;
	}
	
	double x0[14], rvf[6];
	V_Copy(x0, sfpara, 7);
	V_Copy(rvf, &sfpara[7], 6);
	double tf = sfpara[13];
	int outflag = NINT(sfpara[14]);
	// ��ȡMJD_MARSʱ�̻���λ��
	double rv_mars[6];
	V_Copy(rv_mars, &sfpara[15], 6);

	double lam0 = x[0];
	V_Copy(&x0[7], &x[1], 7);
	double chi[4], kappa, dvg[3], tm;
	V_Copy(chi, &x[8], 4);
	kappa = x[12];
	V_Copy(dvg, &x[13], 3);
	tm = x[16];
	



	double dfpara = lam0;

	double AbsTol[14] = {0.0};
	for(int i=0;i<14;i++)
		AbsTol[i]=1e-12;
	double RelTol=1e-12;

	// FILE *fid = NULL;//fopen("temp.txt","w");//����趨��Ч�ļ�·���������Ҫ�ر��ļ�
	// FILE *fid1 = fopen("temp1.txt", "w");
	// FILE *fid2 = fopen("temp2.txt", "w");
	FILE *fid1 = NULL;
	FILE *fid2 = NULL;
	int flag,NumPoint;
	double work[140]={0.0};
	flag = ode45(GA_derivative_type1, x0, &dfpara,  0.0, tm*TOF*86400/TUnit, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid1);
	// ����һ����ƫ��
	// ��������λ��Լ��,7-9
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
	// �����ɳ�����,11
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

	// tm-ʱ���ٶ�Э̬,12-14
	for (int i=0;i<3;++i)
		fvec[12+i] = x0[10+i] - chi[3]*unit_in[i] + kappa*A[i];

	double H1 =  GA_Hamilton_type1(x0, lam0);
	// �����µ�״̬������Э̬����
	for (int i=0;i<3;++i)
	{
		x0[3+i] = x0[3+i] + dvg[i]; // �ٶ�
		x0[7+i] = x0[7+i] - chi[i]; // λ��Э̬
		x0[10+i] = chi[3]*unit_out[i] + kappa*B[i]; // �ٶ�Э̬
	}
	double H2 = GA_Hamilton_type1(x0, lam0);
	// ��̬����ƫ��,15
	double tempu[3];
	V_Minus(tempu, unit_out, unit_in, 3);
	fvec[15] = H1 - H2 - V_Dot(chi, &rv_tm[3], 3) + chi[3]*V_Dot(tempu, am, 3) - kappa*c;


	// �ڶ��׶εĻ���
	// flag = ode45(GA_derivative, x0, dfpara,  tm*TOF*86400/TUnit,  tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);
	flag = ode45(GA_derivative_type1, x0, &dfpara, tm*TOF*86400/TUnit, tf, 14, NumPoint, work, AbsTol, RelTol, 0, -1, -1, fid2);

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

// ��16������²�ֵת��Ϊ17����б���
// [0]:lam0
// [1-7]:7,��ʼʱ��Э̬����
// [8-11]:4,��ʽԼ������
// [12]:����ʽԼ������
// [13-15]:3,�ٶ�����
// [16]:��������ʱ��
void guess2lam(const double* guessValue, double* lamValue)
{
	double temp1, temp2;
	int i;
	// ��1���²�ֵ��Ϊǰ8����б�����norm��Ȼ�����Э̬��һ�����������5�����ӵ�norm
	double norm1, norm2;
	norm1 = guessValue[0];
	norm2 = sqrt(1 - norm1*norm1);
	// ��7���²�ֵ�����norm1������ǰ8����б���
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
	// ��4���²�ֵ�����norm2������5��Լ���ĳ���
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
	// ��1���²�ֵ�������ٶ�������ֵ
	double vAmplitude = sqrt(mpp/rmin)*guessValue[12];
	// ��2���²�ֵ�����vAmplitude������3���ٶȷ���
	double phi, delta;
	phi = (guessValue[13] - 0.5)*DPI;
	delta = guessValue[14]*D2PI;
	lamValue[13] = vAmplitude*cos(phi)*cos(delta);
	lamValue[14] = vAmplitude*cos(phi)*sin(delta);
	lamValue[15] = vAmplitude*sin(phi);
	// ��1���²�ֵ��������������ʱ��
	lamValue[16] = (2.0 + guessValue[15])*365.25/2201;
}

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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa��ά��������544
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
		// ��һ�������ʣ������Ϊ16033.297
		// �����㲻��ʽԼ��
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
		*/

		// �ڶ��������ʣ������Ϊ15360.170
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

		/*
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
		
		// epsi=1.2e-5ʱ�Ľ��
		x[0] = 8.152167e-1;
		x[1] = -2.143450e-1;
		x[2] = -2.975502e-1;
		x[3] = -3.214474e-2;
		x[4] = 2.826663e-1;
		x[5] = -2.109592e-1;
		x[6] = -8.493628e-2;
		x[7] = 1.781236e-1;
		x[8] = 2.541862e-2;
		x[9] = -5.919471e-2;
		x[10] = 7.655084e-2;
		x[11] = -1.616604e-1;
		x[12] = 2.044634e-2;
		x[13] = 6.413660e-2;
		x[14] = 9.066606e-2;
		x[15] = -1.837024e-3;
		x[16] = 3.878463e-1;
		*/

		info = hybrd1(GA_fvec, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()�������־
		// std::cout << enorm(n,fvec) << std::endl;
		if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
		{
			sfpara[15]=1.0;
			int _j = GA_fvec(n, x, fvec, 1, sfpara); // ����ͬ�׼���õ���Э̬��ֵ���������һ�εĻ�����⣨ֱ���������С��ͬ�ײ����������ư����ƵĽ����
			if(fvec[0]>Out[0]) // ʣ������Ϊ�����������о���Ҫ��ֹͣ
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
	printf("���ɹ�%d\n",flag);
	printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
	// printf("lamda0Ϊ:%.6e\n", Out[1]);
	printf("17����б���ֵΪ:\n");
	for (int j=1; j<=17; j++)
		printf("%.6e,\n", Out[j]);
	return flag;
}

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
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa��ά��������544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	// epsi=1.0ʱ�Ľ��
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
	
	// ����һ�������Ϊ��ֵ
	if (Out[0]>1e-8)
		V_Copy(x, &Out[1], 17);

	printf("��ǰͬ�ײ���Ϊ:%f\n", epsi);
	info = hybrd1(GA_fvec, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()�������־
	if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
	{
		sfpara[15]=1.0;
		int _j = GA_fvec(n, x, fvec, 1, sfpara); // ����ͬ�׼���õ���Э̬��ֵ���������һ�εĻ�����⣨ֱ���������С��ͬ�ײ����������ư����ƵĽ����
		if(fvec[0]>Out[0]) // ʣ������Ϊ�����������о���Ҫ��ֹͣ
		{
			flag=1;
			Out[0]=fvec[0];
			V_Copy(&Out[1], x, 17);
		}
		sfpara[15]=0.0;
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
		// printf("lamda0Ϊ:%.6e\n", Out[1]);
		
		printf("17����б���ֵΪ:\n");
		for (int j=1; j<=17; j++)
			printf("%.6e,\n", Out[j]);
		
	}
	else
		printf("���ʧ��%d\n",flag);

	return flag;
}

int GA_FOP_type1(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t)
{
	double sfpara[21] = {0.0};
	V_Copy(sfpara, rv0, 6);
	sfpara[6] = m0;
	V_Copy(&sfpara[7], rv1, 6);
	sfpara[13] = tof;
	sfpara[14] = 0.0;
	V_Copy(&sfpara[15], rv_middle, 6);

	
	int info, flag = 0;
	const int n = 17;
	double x[17] = {0.0}, fvec[17] = {0.0}, wa[600] = {0.0}; // wa��ά��������544
	double guessArray[16] = {0.0};
	double xtol = 1.0e-8;
	
	// epsi=1.2e-5ʱ�Ľ��
	x[0] = 8.152167e-1;
	x[1] = -2.143450e-1;
	x[2] = -2.975502e-1;
	x[3] = -3.214474e-2;
	x[4] = 2.826663e-1;
	x[5] = -2.109592e-1;
	x[6] = -8.493628e-2;
	x[7] = 1.781236e-1;
	x[8] = 2.541862e-2;
	x[9] = -5.919471e-2;
	x[10] = 7.655084e-2;
	x[11] = -1.616604e-1;
	x[12] = 2.044634e-2;
	x[13] = 6.413660e-2;
	x[14] = 9.066606e-2;
	x[15] = -1.837024e-3;
	x[16] = 3.878463e-1;
	

	printf("bang-bang���Ƽ�������\n");
	info = hybrd1(GA_fvec_type1, n, x, fvec, sfpara, wa, xtol, 5, 2000); // info-hybrd1()�������־
	if(info>0 && enorm(n,fvec)<1e-8 && x[0]>0.0)
	{
		sfpara[15]=1.0;
		int _j = GA_fvec_type1(n, x, fvec, 1, sfpara); // ����ͬ�׼���õ���Э̬��ֵ���������һ�εĻ�����⣨ֱ���������С��ͬ�ײ����������ư����ƵĽ����
		if(fvec[0]>Out[0]) // ʣ������Ϊ�����������о���Ҫ��ֹͣ
		{
			flag=1;
			Out[0]=fvec[0];
			V_Copy(&Out[1], x, 17);
		}
		sfpara[14]=0.0;
		printf("���ɹ�%d\n",flag);
		printf("ʣ������Ϊ:%.3fkg\n", Out[0]*MUnit);
		// printf("lamda0Ϊ:%.6e\n", Out[1]);
		
		printf("17����б���ֵΪ:\n");
		for (int j=1; j<=17; j++)
			printf("%.6e,\n", Out[j]);
		
	}
	else
		printf("���ʧ��%d\n",flag);

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
		else if (epsi > 0.00001)
			epsi = epsi - 0.000001;
		else
			break;
		flag = GA_FOP_2(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_middle, PSO_t);
	}
	return flag;
}