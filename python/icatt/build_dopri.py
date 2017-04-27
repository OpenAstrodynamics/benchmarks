from cffi import FFI
import os.path

libpath = os.path.abspath(os.path.join("..","..","lib","dopri","build"))
includepath = os.path.abspath(os.path.join("..","..","lib","dopri"))

ffi = FFI()

ffi.cdef("""
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

if __name__ == "__main__":
    ffi.compile()
