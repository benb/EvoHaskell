#include <math.h>
#include <nlopt.h>
#include <stdio.h>
#include "nlopt_c.h"

int call_nlopt(nlopt_algorithm alg, double xtol, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        nlopt_opt opt;
        double minf; /* the minimum objective value, upon return */
        int retval;
        opt = nlopt_create(alg,np);
        nlopt_set_lower_bounds(opt,lower);
        nlopt_set_upper_bounds(opt,upper);
        nlopt_set_min_objective(opt,func,NULL);
        nlopt_set_xtol_rel(opt,xtol);
        retval = nlopt_optimize(opt, params, &minf) ;
        nlopt_destroy(opt);
        return retval;
}
int opt_bobyqa(double xtol, double *params, unsigned np, double (*func) (unsigned,const double*,double*,void*),double *lower, double *upper){
        return call_nlopt(NLOPT_LN_BOBYQA,xtol,params,np,func,lower,upper);
}

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


