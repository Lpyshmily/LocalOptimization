#include "rk4.h"
#include "Constant.h"

// y:���ֳ�ֵ;t:��ʼʱ��;h:���ֲ���;s:���ֽ��;dim:ά��;
void RK4(void (*Model)(double t, const double* y, double* yp, const double* para), const double* y, double t, double h, double* s, int dim, const double* dfpara)
{
	int i;
	// ���������С�������л���
	if(fabs(h)<=EPSILON)
	{
		for(i=0;i<dim;i++)
			s[i]=y[i];
		return;
	}
	double* k1=new double[dim];
	double* k2=new double[dim];
	Model(t,y,k2, dfpara);
	for(i=0;i<dim;i++)
	{
		s[i]=h*k2[i];
		k1[i]=0.5*h*k2[i]+y[i];
	}
	Model(t+0.5*h,k1,k2, dfpara);
	for(i=0;i<dim;i++)
	{
		s[i]+=2.0*h*k2[i];
		k1[i]=0.5*h*k2[i]+y[i];
	}
	Model(t+0.5*h,k1,k2, dfpara);
	for(i=0;i<dim;i++)
	{
		s[i]+=2.0*h*k2[i];
		k1[i]=h*k2[i]+y[i];
	}
	Model(t+h,k1,k2, dfpara);
	for(i=0;i<dim;i++)
		s[i]=(s[i]+h*k2[i])/6.0+y[i];

	delete[] k1;
	delete[] k2;
}