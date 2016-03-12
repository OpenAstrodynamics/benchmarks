#include <functional>

#ifndef ICATT_KEPLER_H
#define ICATT_KEPLER_H
namespace kepler {
    double period(double sma, double mu);

    double newton(double x0, std::function<double(double)> const &func, std::function<double(double)> const &deriv,
                  int maxiter, double tol);

    double mean2ecc(double M, double ecc);

    double ecc2true(double E, double ecc);
    void benchmark(int times);
}
#endif //ICATT_KEPLER_H
