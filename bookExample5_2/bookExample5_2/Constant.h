#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<cmath>

const double MJD0 = 59534; // 从地球出发的时间 2021-11-16
const double TOF = 2201;
const double MJDf = MJD0 + TOF; // 到达木星时间 2201天之后 61735 2027-11-26
const double MJD_MARS = 60389; // 火星的位置速度对应的时间 2024-3-20

const double Tmax = 2.26;//最大推力，N
const double Isp = 6000.0;//比冲，s
const double MUnit = 20000.0;//探测器质量归一化时的质量单位，kg

const double DPI = 3.1415926535897932384626433832795;//圆周率
const double D2PI = 6.283185307179586476925286766559;//2倍圆周率
const double D2R = 0.017453292519943295769236907684886;//度转弧度时乘以的系数
const double R2D = 57.295779513082320876798154814105;//弧度转度时乘以的系数
const double g0 = 9.80665;//海平面平均重力加速度，m/s^2

const double mu = 1.32712440018e20;//太阳引力系数，m^3/s^2
const double LUnit = 149597870660.0;//归一化时的长度单位，m
const double TUnit = sqrt(LUnit*LUnit*LUnit/mu);//以LUnit为半径的圆轨道周期的1/(2pi)为归一化时的时间单位，s
const double VUnit = LUnit/TUnit;//归一化时的速度单位
const double AUnit = VUnit/TUnit;//归一化时的加速度单位

const double muNU = (mu/LUnit)*(TUnit/LUnit)*(TUnit/LUnit); // 归一化后的引力系数，数值为1
const double Ispg0NU = Isp*g0/VUnit;
const double TmaxNU = Tmax/(MUnit*AUnit);

const double EPSILON = 1.0e-14;//通用的精度要求



// 甩摆行星参数赋值
const double rpp = 3389.9e3/LUnit; // 火星半径
const double rmin = 3889.9e3/LUnit; // 最小甩摆半径
const double mpp = 42828.3e9/mu; // 火星引力系数

const double STEP = 2.0e-3; // 二次同伦同伦到0.0，成功
// const double STEP = 1.0e-4;
#endif