#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include <cmath>

using Eigen::Vector3d;
using std::runtime_error;

namespace lambert {
    struct LambertResults {
        Vector3d v0;
        Vector3d v;
    };

    double c2(double psi) {
        double res;
        auto eps = 1.0;
        if (psi > eps) {
            res = (1 - cos(sqrt(psi))) / psi;
        } else if (psi < -eps) {
            res = (cosh(sqrt(-psi)) - 1) / (-psi);
        } else {
            res = 1.0 / 2.0;
            auto delta = (-psi) / tgamma(2 + 2 + 1);
            auto k = 1;
            while (res + delta != res) {
                res += delta;
                k += 1;
                delta = pow(-psi, k) / tgamma(2*k + 2 + 1);
            }
        }
        return res;
    }

    double c3(double psi) {
        double res;
        auto eps = 1.0;
        if (psi > eps) {
            res = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi));
        } else if (psi < -eps) {
            res = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi));
        } else {
            res = 1.0 / 6.0;
            auto delta = (-psi) / tgamma(2 + 3 + 1);
            int k = 1;
            while (res + delta != res) {
                res += delta;
                k += 1;
                delta = pow(-psi, k) / tgamma(2*k + 3 + 1);
            }
        }
        return res;
    }

    LambertResults lambert(double k, Vector3d r0, Vector3d r, double tof,
            bool shortway=true, int numiter=35, double rtol=1e-8) {
        int t_m;
        if (shortway) {
            t_m = 1;
        } else {
            t_m = -1;
        }
        auto norm_r0 = r0.norm();
        auto norm_r = r.norm();
        auto cos_dnu = r0.dot(r) / (norm_r * norm_r0);

        auto A = t_m * sqrt(norm_r * norm_r0 * (1 + cos_dnu));

        if (A == 0.0) {
            throw runtime_error("Cannot compute orbit, phase angle is 180 degrees");
        }

        auto psi = 0.0;
        auto psi_low = -4*M_PI;
        auto psi_up = 4*M_PI;

        auto count = 0;
        double y;
        while (count < numiter) {
            y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi));
            if (A > 0.0 & y < 0.0) {
                while (y < 0.0) {
                    psi_low = psi;
                    psi = (0.8 * (1.0 / c3(psi)) *
                            (1.0 - (norm_r0 + norm_r) * sqrt(c2(psi)) / A));
                    y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi));
                }
            }
            auto xi = sqrt(y / c2(psi));
            auto tof_new = (pow(xi, 3) * c3(psi) + A * sqrt(y)) / sqrt(k);

            if (fabs((tof_new - tof) / tof) < rtol) {
                break;
            } else {
                count += 1;
                if (tof_new <= tof) {
                    psi_low = psi;
                } else {
                    psi_up = psi;
                }
                psi = (psi_up + psi_low) / 2;
            }
        }
        if (count > numiter) {
            throw runtime_error("Maximum number of iterations reached");
        }
        auto f = 1 - y / norm_r0;
        auto g = A * sqrt(y / k);
        auto gdot = 1 - y / norm_r;
        auto v0 = (r - f * r0) / g;
        auto v = (gdot * r - r0) / g;
        LambertResults results = {v, v0};
        return results;
    }

    void benchmark(int times) {
        auto mu = 3.986004418e5;
        Vector3d r0(5000.0, 10000.0, 2100.0);
        Vector3d r(-14600.0, 2500.0, 7000.0);
        auto tof = 3600.0;

        auto best = std::numeric_limits<double>::infinity();
        auto worst = -std::numeric_limits<double>::infinity();
        double all = 0;
        for (auto i=0; i < times; i++) {
            auto begin = std::chrono::high_resolution_clock::now();

            lambert(mu, r0, r, tof);

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
