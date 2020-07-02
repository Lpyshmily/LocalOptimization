#ifndef GA_FOP_H_
#define GA_FOP_H_

int GA_FOP(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t);
int GA_FOP_2(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t);
int GA_FOP_type1(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t);
int GA_FOP_type2(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t);
int GA_FOP_type2_process(double* Out, const double* rv0, const double* rv1, double m0, double tof, double epsi, int MaxGuessNum, const double* rv_middle, double PSO_t);
int change_epsi(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t);
int change_epsi_type2(double* Out, const double* rv0, const double* rv1, double m0, double tof, int MaxGuessNum, const double* rv_middle, double PSO_t);

#endif