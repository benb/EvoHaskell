#include <math.h>
#include <nlopt.h>
#include <stdio.h>
#include "nlopt_c.h"

int call_nlopt(nlopt_algorithm alg, double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        nlopt_opt opt;
        double minf; /* the minimum objective value, upon return */
        int retval;
        int i;
        opt = nlopt_create(alg,np);
        nlopt_set_lower_bounds(opt,lower);
        nlopt_set_upper_bounds(opt,upper);
        nlopt_set_min_objective(opt,func,NULL);
        nlopt_set_ftol_abs(opt,ftol);
        nlopt_set_initial_step(opt,step_size);
        retval = nlopt_optimize(opt, params, &minf) ;
        nlopt_destroy(opt);
        return retval;
}

int opt_bobyqa(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_BOBYQA,ftol,step_size,params,np,func,lower,upper);
}
int opt_cobyla(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_COBYLA,ftol,step_size,params,np,func,lower,upper);
}
int opt_mma(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_MMA,ftol,step_size,params,np,func,lower,upper);
}
int opt_slsqp(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_SLSQP,ftol,step_size,params,np,func,lower,upper);
}
int opt_newton(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_TNEWTON_PRECOND_RESTART,ftol,step_size,params,np,func,lower,upper);
}
int opt_var1(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_VAR1,ftol,step_size,params,np,func,lower,upper);
}
int opt_var2(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_VAR2,ftol,step_size,params,np,func,lower,upper);
}
int opt_lbfgs(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LD_LBFGS,ftol,step_size,params,np,func,lower,upper);
}
int opt_sbplx(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_SBPLX,ftol,step_size,params,np,func,lower,upper);
}
int opt_neldermead(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_NELDERMEAD,ftol,step_size,params,np,func,lower,upper);
}
int opt_praxis(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_PRAXIS,ftol,step_size,params,np,func,lower,upper);
}
int opt_newuoa(double ftol, double *step_size, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_NEWUOA_BOUND,ftol,step_size,params,np,func,lower,upper);
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


