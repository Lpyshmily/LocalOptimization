#ifndef _MODELFUN_H_
#define _MODELFUN_H_

// ��֪���ֵ�������ee����6���Ƕ�����ֵL0��ĩֵLt�����L0�仯��Lt�����ʱ�䣬musΪ��������ϵ�����������չ������
double L0Lt2dt(double L0, double Lt, const double* ee, double mus);

// ��֪���ֵ�������ee����6���Ƕ�����ֵL0���󾭹�ʱ��dt���Lt��musΪ��������ϵ�����������չ������
double L0dt2Lt(double L0, double dt, const double* ee, double mus);

//�󴺷ֵ�������ee�����¶���ѧ�����е�M����(6*3)��D�����6������ֵD6��
//M����D6��ee��ƫ�������׾���dM(6*3*6)��һ������dD6(6)
//M����D6��ee�Ķ���ƫ�������׾���ddM(6*3*6*6)�Ͷ�������ddD6(6*6)
//step:1,ֻ��M, D6; 2, �㵽dM, dD6;3, �㵽ddM, ddD6
void MatMD(double * M, double& D6, double * dM, double* dD6, double* ddM, double* ddD6, const double* ee, int step);

//ȼ���������⶯��ѧ΢�ַ��̣�����״̬����Э̬������14ά
//tΪʱ��,xΪ��ǰ״̬��Э̬,dxΪ��ʱ��ı仯��, dfpara������Ĳ�����
//����ɹ�ʱ������ֵΪ1�����򷵻�С��1������ֵ
int DynEquFuelOpt(double t, const double* x, double* dx, const double* dfpara);

//����ʱ���״̬�̶���ȼ�����������з���
//sfpara
//[0-6],��ʼ���ֵ���������
//[7-12],ĩ�˴��ֵ���
//[13-14],ת��ʱ���ͬ�ײ���epsilon
//[15],�����Ϣ��ʶ
//x:[0],lam0;[1-7]��ʼЭ̬
//nΪ����ά����fvecΪ�вiflagΪ0ʱ��ʾ�Ϸ����ú�����֪������в�
//����ɹ�ʱ������0�����0��ֵ�����ɹ�ʱ����С��0��ֵ
int ShootFunFuelOpt(int n, const double *x, double *fvec, int iflag, const double* sfpara);

//��������³�ֵ��ͨ��ShootFunFuelOpt()���������ֵ�����µĻ��ֽ������ͨ��hybrid()���д�з������
//�����˹̶���ȼ����������
//OutΪ���������ֱ�Ϊʣ��������8�����ӳ�ֵ����һ��Ϊlamd0)����9ά
//ee0Ϊ�����Ĵ��ֵ���������ee1Ϊ����Ĵ��ֵ���������m0Ϊ����������ftΪ����ʱ�䣬
//epsiΪͬ�ײ�����ֵ��ֱ��ȡһ����С��ֵ����0.001����ýӽ������ƵĽ������MaxGuessNumΪ���²����
//����ɹ�ʱ������1�����򷵻�0
int Fixed2PFuelOpt(double* Out, const double* ee0, const double* ee1, double m0, double ft, double epsi, int MaxGuessNum);

void OutputAmplitude(const double* Out, const double* ee0, double m0, double epsi, const double* TimeNode, int TimeNodeNum);

#endif