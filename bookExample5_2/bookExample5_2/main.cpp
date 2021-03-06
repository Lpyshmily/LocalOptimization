#include <iostream>
#include <ctime>
#include "Constant.h"
#include "GA_PSO.h"
#include "GA_FOP.h"
#include "ExternHeader.h"

int main()
{
	clock_t start, stop;
	start = clock();
	int flag, i;

	// **********
	// 初始条件设定与归一化
	// 初始、末端位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 5.876420e-1, 7.954627e-1, -3.845203e-5, -5.155764, 3.707833, -3.191945e-4 };
	double rv1[6] = { -5.204974, 1.495369, 1.102444e-1, -7.936872e-1, -2.523063, 2.823220e-2 };
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/(365.25*86400)*TUnit;
		rv1[i] = rv1[i]/(365.25*86400)*TUnit;
		// rv0[i] = rv0[i]/D2PI;
		// rv1[i] = rv1[i]/D2PI;
	}
	// 火星经典轨道根数转换为位置速度
	double coe_mars[6] = { 1.52363312, 0.09327933, 1.84785414, 49.48935357, 286.67090811, 328.88755274};
	for (i=2;i<6;++i)
		coe_mars[i]  = coe_mars[i]*D2R;
	double rv_mars[6] = {0.0};
	coe2rv(flag, rv_mars, coe_mars, muNU);

	double m0 = 1.0;
	double tof = 2201*86400/TUnit;
	// **********


	// **********
	// 使用PSO算法搜索最小变轨速度增量对应的引力辅助时刻
	/*
	double para[19];
	V_Copy(para, rv0, 6);
	V_Copy(&para[6], rv1, 6);
	V_Copy(&para[12], rv_mars, 6);
	para[18] = MJD_MARS;
	double x0 = 0.392329395729;
	double dv = GA_PSO_obj(&x0, para);
	printf("脉冲引力辅助速度增量为%.3f\n", dv*VUnit);
	*/

	std::cout << "下面进行PSO搜索" << std::endl;
	double PSO_t = GA_PSO(rv0, rv1, rv_mars, 5);
	printf("引力辅助时刻为%.3f\n", PSO_t);


	
	// 间接法求解引力辅助
	// 求解算法的一些参数设置
	double epsi = 0.0; // 同伦参数设置为1.0
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[18] = {0.0}; // 输出计算结果，0-剩余质量，1~17-17个需要打靶的协态初值
	// flag = GA_FOP(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_mars, PSO_t);
	/*
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	// printf("lamda0为:%.6e\n", Out[1]);
	printf("17个打靶变量值为:\n");
	for (int j=1; j<=17; j++)
		printf("%.6e,\n", Out[j]);
	*/
	// flag = GA_FOP(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_mars, PSO_t); // 对数同伦
	// flag = change_epsi(Out, rv0, rv1, m0, tof, MaxGuessNum, rv_mars, PSO_t); // 对数同伦过程
	flag = GA_FOP_type2(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_mars, PSO_t); // 二次同伦
	flag = GA_FOP_type1(Out, rv0, rv1, m0, tof, MaxGuessNum, rv_mars, PSO_t); // bang-bang控制
	// flag = GA_FOP_type2_process(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum, rv_mars, PSO_t);
	// flag = change_epsi_type2(Out, rv0, rv1, m0, tof, MaxGuessNum, rv_mars, PSO_t); // 二次同伦过程

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);

	return 0;
}