#include <chrono>
#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include "dopri.h"
#include "kepler.h"
#include "elements.h"

using Eigen::VectorXd;
using Eigen::Vector3d;
using std::pow;
using std::sqrt;

namespace dopri {
    void gravity(int *n, double *x, double *y, double *f, double *rpar, int *ipar) {
        auto r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        auto r3 = r*r*r;
        f[0] = y[3];
        f[1] = y[4];
        f[2] = y[5];
        f[3] = -rpar[0] * y[0] / r3;
        f[4] = -rpar[0] * y[1] / r3;
        f[5] = -rpar[0] * y[2] / r3;
    }

    void solout_dummy(int *nr, double *xold, double *x, double *y, int *n, double *con,
                      int *icomp, int *nd, double *rpar, int *ipar, int *irtrn, double *xout){};

    void integrate(void (*func)(int *, double *, double *, double *, double *, int *),
            double *x, VectorXd *rv, double xend, double rpar[], int ipar[],
            double reltol = 1e-6, double abstol = 1e-8) {
        int n = rv->size();
        double rtol[] = {reltol};
        double atol[] = {abstol};
        int itol = 0;
        int iout = 0;
        int lwork = 11*n+8*n+21;
        int liwork = n + 21;
        double work[lwork];
        memset(work, 0, sizeof(work));
        int iwork[liwork];
        memset(iwork, 0, sizeof(iwork));
        int idid = 0;
        c_dop853(&n, func, x, rv->data(), &xend, rtol, atol, &itol, &solout_dummy,
            &iout, work, &lwork, iwork, &liwork, rpar, ipar, &idid);
    }

    void benchmark(int times) {
        auto mu = 3.986004418e5;
        Vector3d r(8.59072560e+02, -4.13720368e+03, 5.29556871e+03);
        Vector3d v(7.37289205e+00, 2.08223573e+00, 4.39999794e-01);
        VectorXd rv(r.size()+v.size());
        rv << r, v;
        VectorXd rv0(rv);
        auto el = elements::elements(r, v, mu);
        auto x = 0.0;
        double rpar[] = {mu};
        int ipar[] = {0};
        auto xend = kepler::period(el[0], mu);
        auto best = std::numeric_limits<double>::infinity();
        auto worst = -std::numeric_limits<double>::infinity();
        double all = 0;
        for (auto i=0; i < times; i++) {
            auto begin = std::chrono::high_resolution_clock::now();
            integrate(&gravity, &x, &rv0, xend, rpar, ipar);
            auto end = std::chrono::high_resolution_clock::now();
            auto current = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e9;
            all += current;
            if (current < best) {
                best = current;
            }
            if (current > worst) {
                worst = current;
            }
            rv0 = rv;
            x = 0;
        }
        std::cout << "[" << all/times << "," << best << "," << worst << "]" << std::endl;
    }
}
