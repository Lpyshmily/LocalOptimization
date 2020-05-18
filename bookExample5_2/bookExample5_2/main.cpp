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
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/D2PI;
		rv1[i] = rv1[i]/D2PI;
	}
	// 火星经典轨道根数转换为位置速度
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
	// 求解Lambert问题，使用PSO算法搜索最小变轨速度增量对应的引力辅助时刻
	double para[19];
	V_Copy(para, rv0, 6);
	V_Copy(&para[6], rv1, 6);
	V_Copy(&para[12], rv_mars, 6);
	para[18] = MJD_MARS;
	double x0 = 0.392329395729;
	double dv = GA_PSO_obj(&x0, para);
	printf("脉冲引力辅助速度增量为%.3f\n", dv*VUnit);

	std::cout << "下面进行PSO搜索" << std::endl;
	dv = GA_PSO(rv0, rv1, MJD_MARS, rv_mars, 5);

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}