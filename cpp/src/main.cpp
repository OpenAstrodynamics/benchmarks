#include "elements.h"
#include "kepler.h"
#include "cppdopri.h"
#include "lambert.h"

int main()
{
    int n = 100000;
    elements::benchmark(n);
    kepler::benchmark(n);
    lambert::benchmark(n);
    dopri::benchmark(n);
}
