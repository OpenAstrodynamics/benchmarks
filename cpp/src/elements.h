#include <Eigen/Dense>

#ifndef ICATT_ELEMENTS_H
#define ICATT_ELEMENTS_H
namespace elements {
    using namespace Eigen;
    VectorXd elements(Vector3d r, Vector3d v, double mu);
    void benchmark(int times);
}
#endif //ICATT_ELEMENTS_H
