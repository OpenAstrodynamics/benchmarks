from cffi import FFI
ffi = FFI()
ffi.cdef("""
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
C = ffi.dlopen("../../../deps/dopri/build/libdopri.dylib")

@ffi.callback("void(int *n, double *x, double *y, double *f, double *rpar, int *ipar)")
def gravity(n, x, y, f, rpar, ipar):
    n = 

