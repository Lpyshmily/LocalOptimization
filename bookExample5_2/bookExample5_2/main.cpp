#include <iostream>
#include <ctime>
#include "Constant.h"
#include "gravity_assist.h"
#include "ExternHeader.h"

int main()
{
	clock_t start, stop;
	start = clock();
	int flag;

	// **********
	// ��ʼ�����趨���һ��
	// ��ʼ��ĩ��λ�ú��ٶȣ���λ�ֱ�ΪAU��AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	// ����ʼ��ĩ�˵��ٶȰ���Constant.h�еĻ�׼���й�һ��
	// ����Constant.h�е�ʱ���һ����׼��һ���һ���Ľ��Ӧ����2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/D2PI;
		rv1[i] = rv1[i]/D2PI;
	}
	// ���Ǿ���������ת��Ϊλ���ٶ�
	double coe_mars[6] = { 1.52363312, 0.09327933, 1.84785414, 49.48935357, 286.67090811, 328.88755274};
	for (int i=2;i<6;++i)
		coe_mars[i]  = coe_mars[i]*D2R;
	double rv_mars[6] = {0.0};
	coe2rv(flag, rv_mars, coe_mars, muNU);

	double m0 = 20000.0;
	m0 = m0/MUnit;
	double tof = 2201*86400/TUnit;
	// **********


	// **********
	// ���Lambert���⣬ʹ��PSO�㷨������С����ٶ�������Ӧ����������ʱ��
	double para[19];
	V_Copy(para, rv0, 6);
	V_Copy(&para[6], rv1, 6);
	V_Copy(&para[12], rv_mars, 6);
	para[18] = MJD_MARS;
	double x0 = 0.392329395729;
	double dv = GA_PSO_obj(&x0, para);
	printf("�������������ٶ�����Ϊ%.3f\n", dv*VUnit);

	std::cout << "�������PSO����" << std::endl;
	dv = GA_PSO(rv0, rv1, MJD_MARS, rv_mars, 5);

	stop = clock();
	printf("������ʱΪ��%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}