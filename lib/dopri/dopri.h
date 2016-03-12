#ifdef __cplusplus
extern "C"
{
#endif

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

void c_dopri5(
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

double c_contd8(int *ii, double *x, double *con, int *icomp, int *nd);
double c_contd5(int *ii, double *x, double *con, int *icomp, int *nd);

#ifdef __cplusplus
}
#endif
