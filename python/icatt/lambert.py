import time
import numpy as np
from math import gamma

from numba import jit, njit

def benchmark(times):
    r0 = np.array([5000.0, 10000.0, 2100.0])
    r = np.array([-14600.0, 2500.0, 7000.0])
    tof = 3600.0
    mu = 3.986004418e5
    lambert(mu, r0, r, tof)

    best = np.inf
    worst = -np.inf
    total = 0.0
    for _ in range(times):
        t0 = time.clock()
        lambert(mu, r0, r, tof)
        t1 = time.clock()
        current = t1 - t0
        if current < best:
            best = current
        if current > worst:
            worst = current
        total += current
    print("[",total/times,",",best,",",worst,"]")

@njit
def lambert(k, r0, r, tof, short=True, numiter=35, rtol=1e-8):
    if short:
        t_m = +1
    else:
        t_m = -1

    norm_r0 = np.dot(r0, r0)**.5
    norm_r = np.dot(r, r)**.5
    cos_dnu = np.dot(r0, r) / (norm_r0 * norm_r)

    A = t_m * (norm_r * norm_r0 * (1 + cos_dnu))**.5

    if A == 0.0:
        raise RuntimeError("Cannot compute orbit, phase angle is 180 degrees")

    psi = 0.0
    psi_low = -4 * np.pi
    psi_up = 4 * np.pi

    count = 0
    while count < numiter:
        y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / c2(psi)**.5
        if A > 0.0 and y < 0.0:
            # Readjust xi_low until y > 0.0
            # Translated directly from Vallado
            while y < 0.0:
                psi_low = psi
                psi = (0.8 * (1.0 / c3(psi)) *
                       (1.0 - (norm_r0 + norm_r) * np.sqrt(c2(psi)) / A))
                y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / c2(psi)**.5

        xi = np.sqrt(y / c2(psi))
        tof_new = (xi**3 * c3(psi) + A * np.sqrt(y)) / np.sqrt(k)

        # Convergence check
        if np.abs((tof_new - tof) / tof) < rtol:
            break
        else:
            count += 1
            # Bisection check
            if tof_new <= tof:
                psi_low = psi
            else:
                psi_up = psi
            psi = (psi_up + psi_low) / 2
    else:
        raise RuntimeError("Maximum number of iterations reached")

    f = 1 - y / norm_r0
    g = A * np.sqrt(y / k)

    gdot = 1 - y / norm_r

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g

    return v0, v

@njit('f8(f8)')
def c2(psi):
    r"""Second Stumpff function.
    For positive arguments:
    .. math::
        c_2(\psi) = \frac{1 - \cos{\sqrt{\psi}}}{\psi}
    """
    eps = 1.0
    if psi > eps:
        res = (1 - np.cos(np.sqrt(psi))) / psi
    elif psi < -eps:
        res = (np.cosh(np.sqrt(-psi)) - 1) / (-psi)
    else:
        res = 1.0 / 2.0
        delta = (-psi) / gamma(2 + 2 + 1)
        k = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / gamma(2 * k + 2 + 1)

    return res


@njit('f8(f8)')
def c3(psi):
    r"""Third Stumpff function.
    For positive arguments:
    .. math::
        c_3(\psi) = \frac{\sqrt{\psi} - \sin{\sqrt{\psi}}}{\sqrt{\psi^3}}
    """
    eps = 1.0
    if psi > eps:
        res = (np.sqrt(psi) - np.sin(np.sqrt(psi))) / (psi * np.sqrt(psi))
    elif psi < -eps:
        res = (np.sinh(np.sqrt(-psi)) - np.sqrt(-psi)) / (-psi * np.sqrt(-psi))
    else:
        res = 1.0 / 6.0
        delta = (-psi) / gamma(2 + 3 + 1)
        k = 1
        while res + delta != res:
            res = res + delta
            k += 1
            delta = (-psi) ** k / gamma(2 * k + 3 + 1)

    return res
