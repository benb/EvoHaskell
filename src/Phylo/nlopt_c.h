int opt_bobyqa(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper);
int opt_cobyla(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper);
int opt_mma(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper);
