#ifndef _GRAVITY_ASSIST_H_
#define _GRAVITY_ASSIST_H_

double GA_PSO_obj(const double* px, const double* para);
double GA_PSO(const double* rv0, const double* rv1, double mjd_middle, const double* rv_middle, int psotime);

#endif