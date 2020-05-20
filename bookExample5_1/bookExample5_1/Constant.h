#ifndef _CONSTANT_H_
#define _CONSTANT_H_
#include<math.h>

// ���������ֵ���ȳ塢��ʼ����
// ������Ҫ���йأ�ÿ����Ҫ�޸�
const double Tmax=0.33;//���������N
const double Isp=3800.0;//�ȳ壬s
const double MUnit=1000.0;//̽����������һ��ʱ��������λ��kg

// �����һ����׼��һЩ����������Ҫ�޸�
const double mu=1.32712440018e20;//̫������ϵ����m^3/s^2
const double LUnit=149597870660.0;//��һ��ʱ�ĳ��ȵ�λ��m
const double TUnit=sqrt(LUnit*LUnit*LUnit/mu);//��LUnitΪ�뾶��Բ������ڵ�1/(2pi)Ϊ��һ��ʱ��ʱ�䵥λ��s
const double VUnit=LUnit/TUnit;//��һ��ʱ���ٶȵ�λ
const double AUnit=VUnit/TUnit;//��һ��ʱ�ļ��ٶȵ�λ

const double DPI=3.1415926535897932384626433832795;//Բ����
const double D2PI=6.283185307179586476925286766559;//2��Բ����
const double D2R=0.017453292519943295769236907684886;//��ת����ʱ���Ե�ϵ��
const double R2D=57.295779513082320876798154814105;//����ת��ʱ���Ե�ϵ��
const double g0=9.80665;//��ƽ��ƽ���������ٶȣ�m/s^2

const double muNU=(mu/LUnit)*(TUnit/LUnit)*(TUnit/LUnit); // ��һ�����̫������ϵ����ֵΪ1
const double Ispg0NU=Isp*g0/VUnit;
const double TmaxNU=Tmax/(MUnit*AUnit);

const double EPSILON=1.0e-14;//ͨ�õľ���Ҫ��
#endif