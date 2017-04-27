import numpy as np
from numba import jit, njit
import time
from icatt import elements

def benchmark(times):
    r = np.array([8.59072560e+02, -4.13720368e+03, 5.29556871e+03])
    v = np.array([7.37289205e+00, 2.08223573e+00, 4.39999794e-01])
    tof = 3600.0
    mu = 3.986004418e5
    el = elements.elements(r, v, mu)
    mean2ecc(np.pi/2, el[1])

    best = np.inf
    worst = -np.inf
    total = 0.0
    for _ in range(times):
        t0 = time.clock()
        mean2ecc(np.pi/2, el[1])
        t1 = time.clock()
        current = t1 - t0
        if current < best:
            best = current
        if current > worst:
            worst = current
        total += current
    print("[",total/times,",",best,",",worst,"]")

# @jit
def newton(x0, func, derivative, maxiter=50, tol=1e-8):
    p0 = x0
    for _ in range(maxiter):
        p = p0 - func(p0)/derivative(p0)
        if np.abs(p - p0) < tol:
            return p
        p0 = p
    raise RuntimeError("Not converged.")

# @jit
def mean2ecc(M, ecc):
    def keplereq(E):
        return E - ecc*np.sin(E) - M
    def keplerderiv(E):
        return 1 - ecc*np.cos(E)
    return newton(M, keplereq, keplerderiv)

@njit
def period(a, mu):
    return 2*np.pi*np.sqrt(a**3/mu)

