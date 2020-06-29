#ifndef RK4_H_
#define RK4_H_

void RK4(void (*Model)(double t, const double* y, double* yp, const double* para), const double* y, double t, double h, double* s, int dim, const double* dfpara);

#endif