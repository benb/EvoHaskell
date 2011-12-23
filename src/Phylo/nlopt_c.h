int opt_bobyqa(double xtol, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper);
int opt_cobyla(double xtol, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper);
