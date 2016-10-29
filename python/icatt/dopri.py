import sys
import numpy as np
from icatt._dopri import ffi, lib
import icatt.kepler as kepler
import icatt.elements as elements
import time

def dopri_benchmark(times):
    y = np.array([8.59072560e+02, -4.13720368e+03, 5.29556871e+03, 7.37289205e+00, 2.08223573e+00, 4.39999794e-01])
    y0 = y.copy()
    x = 0.0
    mu = 3.986004418e5
    rpar = np.array([mu])
    el = elements.elements(y[0:3], y[3:], mu)
    tp = kepler.period(el[0], mu)
    # dopri(lib.gravity, x, y, tp, rpar)
    dopri(ffi.addressof(lib, "c_gravity"), x, y, tp, rpar)

    best = np.inf
    worst = -np.inf
    total = 0.0
    for _ in range(times):
        t0 = time.clock()
        # dopri(lib.gravity, x, y, tp, rpar)
        dopri(ffi.addressof(lib, "c_gravity"), x, y, tp, rpar)
        t1 = time.clock()
        x = 0.0
        y = y0.copy()
        current = t1 - t0
        if current < best:
            best = current
        if current > worst:
            worst = current
        total += current
    print("[",total/times,",",best,",",worst,"]")

@ffi.def_extern()
def gravity(n, x, y, f, rpar, ipar):
    r = np.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])
    r3 = r*r*r
    f[0] = y[3]
    f[1] = y[4]
    f[2] = y[5]
    f[3] = -rpar[0]*y[0]/r3
    f[4] = -rpar[0]*y[1]/r3
    f[5] = -rpar[0]*y[2]/r3

@ffi.def_extern()
def solout(nr, xold, x, y, n, con, icomp, nd, rpar, ipar, irtrn, xout):
    pass

def dopri(func, x, y, xend, rpar, reltol=1e-6, abstol=1e-8):
    n = len(y)
    lwork = 11*n + 8*n + 21
    liwork = n + 21
    lwork_ = ffi.new("int *", lwork)
    liwork_ = ffi.new("int *", liwork)
    work = np.zeros(lwork)
    iwork = np.zeros(liwork, dtype=np.int32)
    rtol = np.array([reltol])
    atol = np.array([abstol])

    _n  = ffi.new("int *", n)
    _x = ffi.new('double *', x)
    _xend = ffi.new('double *', xend)
    _iout = ffi.new("int *", 0)
    _idid = ffi.new("int *", 0)
    _itol = ffi.new("int *", 0)
    _lwork = ffi.new("int *", lwork)
    _liwork = ffi.new("int *", liwork)
    _ipar = ffi.new("int []", [])
    _y = ffi.cast('double *', y.ctypes.data)
    _work = ffi.cast('double *', work.ctypes.data)
    _iwork = ffi.cast('int *', iwork.ctypes.data)
    _rtol = ffi.cast('double *', rtol.ctypes.data)
    _atol = ffi.cast('double *', atol.ctypes.data)
    _rpar = ffi.cast('double *', rpar.ctypes.data)
    lib.c_dop853(_n, func, _x, _y, _xend, _rtol, _atol, _itol, lib.solout, _iout, _work, _lwork, _iwork, _liwork, _rpar, _ipar, _idid)
