#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "elements.h"

using Eigen::Vector3d;

namespace elements
{
    VectorXd elements(Vector3d r, Vector3d v, double mu)
    {
        auto r_mag = r.norm();
        auto v_mag = v.norm();
        auto v_mag2 = v_mag * v_mag;
        auto h = r.cross(v);
        auto h_mag = h.norm();
        Vector3d k(0,0,1);
        auto n = k.cross(h);
        auto n_mag = n.norm();
        auto xi = v_mag2/2 - mu/r_mag;
        auto e = ((v_mag2-mu/r_mag)*r - v*r.dot(v))/mu;
        auto ecc = e.norm();
        double sma;
        if (ecc != 1) {
            sma = -mu/(2*xi);
        } else {
            sma = pow(h_mag,2)/mu;
        }
        auto inc = acos(h.z()/h_mag);
        auto node = acos(n.x()/n_mag);
        auto peri = acos(n.dot(e)/(ecc*n_mag));
        auto ano = acos(e.dot(r)/(ecc*r_mag));
        if (n.y() < 0) {
            node = M_PI*2 - node;
        }
        if (e.z() < 0) {
            peri = M_PI*2 - peri;
        }
        if (r.dot(v) < 0) {
            ano = M_PI*2 - ano;
        }
        VectorXd out(6); out << sma, ecc, inc, node, peri, ano;
        return out;
    }

    void benchmark(int times) {
        auto mu = 3.986004418e5;
        Vector3d r(8.59072560e+02, -4.13720368e+03, 5.29556871e+03);
        Vector3d v(7.37289205e+00, 2.08223573e+00, 4.39999794e-01);

        auto best = std::numeric_limits<double>::infinity();
        auto worst = -std::numeric_limits<double>::infinity();
        double all = 0;
        for (auto i=0; i < times; i++) {
            auto begin = std::chrono::high_resolution_clock::now();

            elements(r, v, mu);

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
