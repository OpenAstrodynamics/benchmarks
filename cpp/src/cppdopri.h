#include <Eigen/Dense>
#include "dopri.h"

#ifndef ICATT_DOPRI_H
#define ICATT_DOPRI_H
namespace dopri {
    using Eigen::VectorXd;
    void gravity(int *n, double *x, double *y, double *f, double *rpar, int *ipar);

    void solout_dummy(int *nr, double *xold, double *x, double *y, int *n, double *con,
                      int *icomp, int *nd, double *rpar, int *ipar, int *irtrn, double *xout);
    void integrate(void (*func)(int *, double *, double *, double *, double *, int *),
            double *x, VectorXd *rv, double xend, double rpar[], int ipar[],
            double reltol, double abstol);
    void benchmark(int times);
}
#endif //ICATT_DOPRI_H
