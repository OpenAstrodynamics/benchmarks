from cffi import FFI
import os.path

libpath = os.path.abspath(os.path.join("..","..","lib","dopri","build"))
includepath = os.path.abspath(os.path.join("..","..","lib","dopri"))

ffi = FFI()

ffi.cdef("""
extern "Python" void gravity(int *, double *, double *, double *, double *, int *);

void c_gravity(int *n, double *x, double *y, double *f, double *rpar, int *ipar);

extern "Python" void solout(int *, double *, double *, double *, int *, double *,
                      int *, int *, double *, int *, int *, double *);

void c_dop853(
    int *n,
    void (*fcn)(int *n, double *x, double *y, double *f, double *rpar, int *ipar),
    double *x,
    double *y,
    double *xend,
    double *rtol,
    double *atol,
    int *itol,
    void (*solout)(int *nr, double *xold, double *x, double *y, int *n, double *con,
        int *icomp, int *nd, double *rpar, int *ipar, int *irtrn, double *xout),
    int *iout,
    double *work,
    int *lwork,
    int *iwork,
    int *liwork,
    double *rpar,
    int *ipar,
    int *idid);
""")

ffi.set_source("_dopri",
    """
    #include <math.h>
    #include "dopri.h"

    void c_gravity(int *n, double *x, double *y, double *f, double *rpar, int *ipar) {
        double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        double r3 = r*r*r;
        f[0] = y[3];
        f[1] = y[4];
        f[2] = y[5];
        f[3] = -rpar[0] * y[0] / r3;
        f[4] = -rpar[0] * y[1] / r3;
        f[5] = -rpar[0] * y[2] / r3;
    }
    """,
    libraries=["dopri"],
    include_dirs=[includepath],
    library_dirs=[libpath],
    extra_link_args=['-Wl,-rpath,{libpath}'.format(libpath=libpath)],
)

if __name__ == "__main__":
    ffi.compile()
