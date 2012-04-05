#include <math.h>
#include <nlopt.h>
#include <stdio.h>
#include "nlopt_c.h"

int call_nlopt(nlopt_algorithm alg, double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        nlopt_opt opt;
        double minf; /* the minimum objective value, upon return */
        int retval;
        int i;
        opt = nlopt_create(alg,np);
        nlopt_set_lower_bounds(opt,lower);
        nlopt_set_upper_bounds(opt,upper);
        nlopt_set_min_objective(opt,func,NULL);
        nlopt_set_xtol_abs1(opt,xtol);
        nlopt_set_initial_step(opt,step_size);
        retval = nlopt_optimize(opt, params, &minf) ;
        nlopt_destroy(opt);
        return retval;
}

int opt_bobyqa(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_BOBYQA,xtol,step_size,params,np,func,lower,upper);
}
int opt_cobyla(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_COBYLA,xtol,step_size,params,np,func,lower,upper);
}
int opt_mma(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_MMA,xtol,step_size,params,np,func,lower,upper);
}
int opt_slsqp(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_SLSQP,xtol,step_size,params,np,func,lower,upper);
}
int opt_newton(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_TNEWTON_PRECOND_RESTART,xtol,step_size,params,np,func,lower,upper);
}
int opt_var1(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_VAR1,xtol,step_size,params,np,func,lower,upper);
}
int opt_var2(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_VAR2,xtol,step_size,params,np,func,lower,upper);
}
int opt_lbfgs(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_LBFGS,xtol,step_size,params,np,func,lower,upper);
}
int opt_sbplx(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_SBPLX,xtol,step_size,params,np,func,lower,upper);
}
int opt_neldermead(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_NELDERMEAD,xtol,step_size,params,np,func,lower,upper);
}
int opt_praxis(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_PRAXIS,xtol,step_size,params,np,func,lower,upper);
}
int opt_newuoa(double xtol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_NEWUOA_BOUND,xtol,step_size,params,np,func,lower,upper);
}
/*
 * NLOPT_LD_TNEWTON_PRECOND_RESTART
 * NLOPT_LD_VAR2
 * NLOPT_LD_VAR1
 * NLOPT_LD_LBFGS
 * NLOPT_LN_SBPLX
 * NLOPT_LN_NELDERMEAD
 * NLOPT_LN_PRAXIS
 * NLOPT_LN_NEWUOA_BOUND
 *
/*int main(int argc, char* argv[]){*/
        /*double lb[2] = { -5.12, -5 }; [> lower bounds <]*/
        /*double ub[2] = { 5, 5};*/
        /*double x[2] = { 1.234, 3.4 };  [> some initial guess <]*/
        /*double minf; [> the minimum objective value, upon return <]*/
        /*int retval = call_nlopt(NLOPT_LN_BOBYQA,1E-04,x,2,myfunc,lb,ub);*/
        /*printf("ANS %f %f -> %f (%d)\n",x[0],x[1],myfunc(2,x,NULL,NULL),retval);*/
/*}*/


/*double myfunc(unsigned n, const double *x, double *grad, void *my_func_data) {*/
        /*double ans;*/
        /*ans = x[0]*x[0] + (x[1]*(x[1]+1));*/
        /*printf("%f,%f,%f \n",x[0],x[1],ans);*/
        /*return ans;*/
/*}*/


