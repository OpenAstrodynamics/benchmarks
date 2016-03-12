import numpy as np
from numba import njit, jit
import time

def kepler_elements(times):
    r = np.array([8.59072560e+02, -4.13720368e+03, 5.29556871e+03])
    v = np.array([7.37289205e+00, 2.08223573e+00, 4.39999794e-01])
    mu = 3.986004418e5
    el = elements(r, v, mu)

    best = np.inf
    worst = -np.inf
    total = 0.0
    for _ in range(times):
        t0 = time.clock()
        el = elements(r, v, mu)
        t1 = time.clock()
        current = t1 - t0
        if current < best:
            best = current
        if current > worst:
            worst = current
        total += current
    print("[",total/times,",",best,",",worst,"]")

@njit
def elements(r, v, mu):
    k = np.zeros(3)
    k[2] = 1.0
    el = np.empty(6)
    r_mag = np.dot(r, r)**(1/2)
    v_mag = np.dot(v, v)**(1/2)
    h = cross(r, v)
    h_mag = np.dot(h, h)**(1/2)
    n = cross(k, h)
    n_mag = np.dot(n, n)**(1/2)
    xi = v_mag ** 2 / 2 - mu / r_mag
    e = ((v_mag ** 2 - mu / r_mag) * r - v * np.dot(r, v)) / mu
    el[1] = np.dot(e, e)**(1/2)
    if not el[1] == 1:
        el[0] = - mu / (2 * xi)
    else:
        el[0] = h_mag ** 2 / mu
    el[2] = np.arccos(h[2] / h_mag)
    el[3] = np.arctan2(n[1]/h_mag, n[0]/h_mag)
    el[4] = np.arccos(np.dot(n, e) / (el[1] * n_mag))
    el[5] = np.arccos(np.dot(e, r) / (el[1] * r_mag))
    # Quadrant checks
    if n[1] < 0:
        el[3] = 2*np.pi - el[3]
    if e[2] < 0:
        el[4] = 2*np.pi - el[4]
    if np.dot(r, v) < 0:
        el[5] = 2*np.pi - el[5]
    return el

@njit
def cross(a, b):
    c = np.zeros(3)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c
