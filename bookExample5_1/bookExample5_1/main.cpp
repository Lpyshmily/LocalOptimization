#include <iostream>
#include <ctime>
#include "Constant.h"
#include "ExternHeader.h"
#include "ModelFun.h"
#include "Cartesian.h"

int main()
{
	clock_t start, stop;
	start = clock();

	double m0 = 1000.0; // 初始质量，单位kg
	m0 = m0/MUnit; // 质量归一化
	double tof = 1000.0*86400.0/TUnit; // 飞行时间为1000天，进行归一化
	// 初始位置和速度，单位分别为AU和AU/a
	double rv0[6] = { 9.708322e-1, 2.375844e-1,	-1.671055e-6, -1.598191, 6.081958, 9.443368e-5 };
	// 末端位置和速度,，单位分别为AU和AU/a
	double rv1[6] = { -3.277178e-1,	6.389172e-1, 2.765929e-2, -6.598211, -3.412933,	3.340902e-1	};
	// 将初始和末端的速度按照Constant.h中的基准进行归一化
	// 根据Constant.h中的时间归一化基准，一年归一化的结果应该是2pi
	for (int i=3;i<6;++i)
	{
		rv0[i] = rv0[i]/D2PI;
		rv1[i] = rv1[i]/D2PI;
	}
	// 至此，所有初始条件归一化完毕

	
	// 先用春分点轨道根数进行求解
	int flag;
	// 将位置速度转化为春分点轨道根数
	double ee0[6]={0.0}, ee1[6]={0.0};
	rv2ee(flag, ee0, rv0, muNU);
	rv2ee(flag, ee1, rv1, muNU);
	ee1[5] = ee1[5] + 3*D2PI;
	// 求解算法的一些参数设置
	double epsi = 1.0E-5;//取定一个较小的同伦参数直接求解近似邦邦控制的结果
	int MaxGuessNum = 100;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[9] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = Fixed2PFuelOpt(Out, ee0, ee1, m0, tof, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0为:%.6e\n", Out[1]);
	printf("7个初始协态变量值为:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
	

	/*
	// 用直角坐标系进行求解
	int flag;
	// 求解算法的一些参数设置
	double epsi = 1.0E-5;//取定一个较小的同伦参数直接求解近似邦邦控制的结果
	int MaxGuessNum = 500;//设置最大随机猜测次数
	srand( (unsigned)time( NULL ) );//设定随机数种子，若没有此设置，每次产生一样的随机数
	// 求解
	double Out[9] = {0.0}; // 输出计算结果，0-剩余质量，1~8-8个需要打靶的协态初值
	flag = FOP(Out, rv0, rv1, m0, tof, epsi, MaxGuessNum);
	printf("求解成功%d\n",flag);
	printf("剩余质量为:%.3fkg\n", Out[0]*MUnit);
	printf("lamda0为:%.6e\n", Out[1]);
	printf("7个初始协态变量值为:\n");
	for (int j=2; j<=8; j++)
		printf("%.6e,\n", Out[j]);
	*/

	// 离散积分输出推力大小随时间的变化
	if (flag>0)
	{
		int TimeNodeNum = 1000;
		double* TimeNode = new double[TimeNodeNum];
		double seg = tof/TimeNodeNum;
		for (int i=0;i<TimeNodeNum;++i)
			TimeNode[i] = seg*i;
		TimeNode[TimeNodeNum-1] = tof;
		OutputAmplitude(Out, ee0, m0, epsi, TimeNode, TimeNodeNum);
		delete[] TimeNode;
	}

	stop = clock();
	printf("计算用时为：%.3fs\n", (double)(stop-start)/CLOCKS_PER_SEC);
	return 0;
}