#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<cmath>

const double MJD0 = 59534; // �ӵ��������ʱ�� 2021-11-16
const double TOF = 2201;
const double MJDf = MJD0 + TOF; // ����ľ��ʱ�� 2201��֮�� 61735 2027-11-26
const double MJD_MARS = 60389; // ���ǵ�λ���ٶȶ�Ӧ��ʱ�� 2024-3-20

const double Tmax = 2.26;//���������N
const double Isp = 6000.0;//�ȳ壬s
const double MUnit = 20000.0;//̽����������һ��ʱ��������λ��kg

const double DPI = 3.1415926535897932384626433832795;//Բ����
const double D2PI = 6.283185307179586476925286766559;//2��Բ����
const double D2R = 0.017453292519943295769236907684886;//��ת����ʱ���Ե�ϵ��
const double R2D = 57.295779513082320876798154814105;//����ת��ʱ���Ե�ϵ��
const double g0 = 9.80665;//��ƽ��ƽ���������ٶȣ�m/s^2

const double mu = 1.32712440018e20;//̫������ϵ����m^3/s^2
const double LUnit = 149597870660.0;//��һ��ʱ�ĳ��ȵ�λ��m
const double TUnit = sqrt(LUnit*LUnit*LUnit/mu);//��LUnitΪ�뾶��Բ������ڵ�1/(2pi)Ϊ��һ��ʱ��ʱ�䵥λ��s
const double VUnit = LUnit/TUnit;//��һ��ʱ���ٶȵ�λ
const double AUnit = VUnit/TUnit;//��һ��ʱ�ļ��ٶȵ�λ

const double muNU = (mu/LUnit)*(TUnit/LUnit)*(TUnit/LUnit); // ��һ���������ϵ������ֵΪ1
const double Ispg0NU = Isp*g0/VUnit;
const double TmaxNU = Tmax/(MUnit*AUnit);

const double EPSILON = 1.0e-14;//ͨ�õľ���Ҫ��



// ˦�����ǲ�����ֵ
const double rpp = 3389.9e3/LUnit; // ���ǰ뾶
const double rmin = 3889.9e3/LUnit; // ��С˦�ڰ뾶
const double mpp = 42828.3e9/mu; // ��������ϵ��

const double STEP = 2.0e-3; // ����ͬ��ͬ�׵�0.0���ɹ�
// const double STEP = 1.0e-4;
#endif