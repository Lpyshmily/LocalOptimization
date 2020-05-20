#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<math.h>

// 最大推力幅值、比冲、初始质量
// 与任务要求有关，每次需要修改
const double Tmax=0.33;//最大推力，N
const double Isp=3800.0;//比冲，s
const double MUnit=1000.0;//探测器质量归一化时的质量单位，kg

// 定义归一化基准和一些常数，不需要修改
const double mu=1.32712440018e20;//太阳引力系数，m^3/s^2
const double LUnit=149597870660.0;//归一化时的长度单位，m
const double TUnit=sqrt(LUnit*LUnit*LUnit/mu);//以LUnit为半径的圆轨道周期的1/(2pi)为归一化时的时间单位，s
const double VUnit=LUnit/TUnit;//归一化时的速度单位
const double AUnit=VUnit/TUnit;//归一化时的加速度单位

const double DPI=3.1415926535897932384626433832795;//圆周率
const double D2PI=6.283185307179586476925286766559;//2倍圆周率
const double D2R=0.017453292519943295769236907684886;//度转弧度时乘以的系数
const double R2D=57.295779513082320876798154814105;//弧度转度时乘以的系数
const double g0=9.80665;//海平面平均重力加速度，m/s^2

const double muNU=(mu/LUnit)*(TUnit/LUnit)*(TUnit/LUnit); // 归一化后的太阳引力系数，值为1
const double Ispg0NU=Isp*g0/VUnit;
const double TmaxNU=Tmax/(MUnit*AUnit);

const double EPSILON=1.0e-14;//通用的精度要求
#endif