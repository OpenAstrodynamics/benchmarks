#include <math.h>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>
#include "kepler.h"
#include "elements.h"

using std::function;
using std::runtime_error;
using Eigen::Vector3d;

namespace kepler {
    double period(double sma, double mu) {
        return M_PI * 2 * sqrt(pow(sma, 3) / mu);
    }

    double newton(
            double p0,
            function<double(double)> const &func,
            function<double(double)> const &deriv,
            int maxiter = 50,
            double tol = 1e-8
    ) {
        for (auto i = 1; i < maxiter; i++) {
            auto p = p0 - func(p0) / deriv(p0);
            if (fabs(p - p0) < tol) {
                return p;
            }
            p0 = p;
        }
        throw runtime_error("Not converged.");
    }

    double mean2ecc(double M, double ecc) {
        auto E = newton(M, [ecc, M](double E) -> double {
            return E - ecc * sin(E) - M;
        }, [ecc](double E) -> double {
            return 1 - ecc * cos(E);
        });
        return E;
    }

    double ecc2true(double E, double ecc) {
        return 2 * atan2(sqrt(1 + ecc) * sin(E / 2), sqrt(1 - ecc) * cos(E / 2));
    }

    void benchmark(int times) {
        auto mu = 3.986004418e5;
        Vector3d r(8.59072560e+02, -4.13720368e+03, 5.29556871e+03);
        Vector3d v(7.37289205e+00, 2.08223573e+00, 4.39999794e-01);
        auto el = elements::elements(r, v, mu);

        auto best = std::numeric_limits<double>::infinity();
        auto worst = -std::numeric_limits<double>::infinity();
        double all = 0;
        for (auto i=0; i < times; i++) {
            auto begin = std::chrono::high_resolution_clock::now();

            mean2ecc(M_PI, el[1]);

            auto end = std::chrono::high_resolution_clock::now();
            auto current = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e9;
            all += current;
            if (current < best) {
                best = current;
            }
            if (current > worst) {
                worst = current;
            }
        }
        std::cout << "[" << all/times << "," << best << "," << worst << "]" << std::endl;
    }
}

