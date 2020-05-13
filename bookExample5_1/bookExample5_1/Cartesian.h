#ifndef _CARTESIAN_H_
#define _CARTESIAN_H_

int Derivative(double t, const double* x, double* dx, const double* dfpara);
int ComputeFvec(int n, const double *x, double *fvec, int iflag, const double* sfpara);
int FOP(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum);

#endif