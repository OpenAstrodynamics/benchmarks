/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include <string.h>
#include "mex.h"
#include "dopri.h"

char* g_buf;
void solout_dummy(int *nr, double *xold, double *x, double *y, int *n, double *con,
      int *icomp, int *nd, double *rpar, int *ipar, int *irtrn, double *xout) { };

void gravity(int *n, double *x, double *y, double *f, double *rpar, int *ipar) {
    int nlhs = 2;
    int nrhs = 1;
    double *ypointer;
    double *mupointer;
    double *fpointer;
    mxArray *plhs[nlhs];
    mxArray *prhs[nrhs];
    plhs[0] = mxCreateDoubleMatrix(*n, 0, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 0, mxREAL);
    prhs[0] = mxCreateDoubleMatrix(*n, 0, mxREAL);
    ypointer = mxGetPr(plhs[0]);
    mupointer = mxGetPr(plhs[1]);
    fpointer = mxGetPr(prhs[0]);
    mupointer[0] = rpar[0];
    for (int i=0; i < 6; i++) {
        ypointer[i] = y[i];
    }
    mexErrMsgIdAndTxt("icatt:propagator:nlhs","One output required.");
    int ret = mexCallMATLAB(nlhs, plhs, nrhs, prhs, g_buf);
    if (ret != 0) {
        mexErrMsgIdAndTxt("icatt:propagator:gravity","Callback failed.");
    }
    for (int i=0; i < 6; i++) {
        f[i] = fpointer[i];
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char* buf;
    double x;
    double *y;
    /* size_t n; */
    double *y1;
    double xend;
    double *rpar;
    double *atol;
    double *rtol;

    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("icatt:propagator:nrhs","Seven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("icatt:propagator:nlhs","One output required.");
    }
    /* make sure the first input argument is a string */
    if( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notChar","First argument must be a string.");
    }
    /* make sure the second input argument is a scalar */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }

    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }

    mwSize buflen;
    buflen = mxGetNumberOfElements(prhs[0]) + 1;
    g_buf = mxCalloc(buflen, sizeof(char));
    if (mxGetString(prhs[0], g_buf, buflen) != 0)
        mexErrMsgIdAndTxt( "MATLAB:explore:invalidStringArray",
                "Could not convert string data.");
    /* get the value of the scalar input  */
    x = mxGetScalar(prhs[1]);

    /* create a pointer to the real data in the input matrix  */
    y = mxGetPr(prhs[2]);

    /* get dimensions of the input matrix */
    /* n = mxGetN(prhs[3]); */
    /* int m = (int) n; */

    int n = 6;
    xend = mxGetScalar(prhs[3]);
    rtol = mxGetPr(prhs[4]);
    atol = mxGetPr(prhs[5]);
    rpar = mxGetPr(prhs[6]);
    int itol = 0;
    int iout = 0;
    int lwork = 11*n+8*n+21;
    int liwork = n + 21;
    double *work;
    int *iwork;
    int *ipar;
    work = (double*) mxMalloc(lwork*sizeof(double));
    iwork = (int*) mxMalloc(liwork*sizeof(int));
    ipar = (int*) mxMalloc(sizeof(int));
    /* double work[lwork]; */
    /* memset(work, 0, sizeof(work)); */
    /* int iwork[liwork]; */
    /* memset(iwork, 0, sizeof(iwork)); */
    int idid = 0;
    /* int ipar[] = {0}; */
    /*  */
    c_dop853(&n, &gravity, &x, y, &xend, rtol, atol, &itol, &solout_dummy,
        &iout, work, &lwork, iwork, &liwork, rpar, ipar, &idid);


    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)n,mxREAL);

    /* get a pointer to the real data in the output matrix */
    y1 = mxGetPr(plhs[0]);

}

