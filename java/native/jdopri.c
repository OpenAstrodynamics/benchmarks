#include <jni.h>
#include <stdio.h>
#include <string.h>
#include "dopri.h"

jobject g_DopriInterface;
jmethodID g_FuncMethod;
JNIEnv *g_env;

void fcn(int *n, double *x, double *y, double *f, double *rpar, int *ipar){
    int l = ipar[0];
    jdoubleArray jrpar = (*g_env)->NewDoubleArray(g_env, l);
    (*g_env)->SetDoubleArrayRegion(g_env, jrpar, 0 , l, rpar);
    jdoubleArray jy = (*g_env)->NewDoubleArray(g_env, *n);
    (*g_env)->SetDoubleArrayRegion(g_env, jy, 0 , *n, y);

    jdoubleArray retval;
    retval = (*g_env)->CallObjectMethod(g_env, g_DopriInterface, g_FuncMethod, *x, jy, jrpar);
    jdouble *out = (*g_env)->GetDoubleArrayElements(g_env, retval, NULL);
    memcpy(f, out, *n*sizeof(double));
    (*g_env)->ReleaseDoubleArrayElements(g_env, retval, out, 0);
    return;
}

void solout(int *nr, double *xold, double *x, double *y, int *n, double *con,
        int *icomp, int *nd, double *rpar, int *ipar, int *irtrn, double *xout){
    return;
}

JNIEXPORT void JNICALL Java_com_helgeeichhorn_icatt_Dopri_jdop853
  (JNIEnv *env, jobject thisObj, jobject DopriInterface, jint n, jdouble x, jdoubleArray y,
   jdouble xend, jdoubleArray rtol, jdoubleArray atol, jint itol, jint iout,
   jdoubleArray rpar, jintArray ipar, jint idid) {
    jclass objClass = (*env)->GetObjectClass(env, DopriInterface);
    jmethodID FuncMethod = (*env)->GetMethodID(env, objClass, "func", "(D[D[D)[D");
    if (FuncMethod == NULL) {
        printf("'func' method not found.\n");
        return;
    }
    g_FuncMethod = FuncMethod;
    g_DopriInterface = DopriInterface;
    g_env = env;
    jdouble *cy = (*env)->GetDoubleArrayElements(env, y, NULL);
    jdouble *crpar = (*env)->GetDoubleArrayElements(env, rpar, NULL);
    jdouble *crtol = (*env)->GetDoubleArrayElements(env, rtol, NULL);
    jdouble *catol = (*env)->GetDoubleArrayElements(env, atol, NULL);
    jint *cipar = (*env)->GetIntArrayElements(env, ipar, NULL);
    int lwork = 11*n+8*n+21;
    int liwork = n + 21;
    double work[lwork];
    memset(work, 0, sizeof(work));
    int iwork[liwork];
    memset(iwork, 0, sizeof(iwork));
    c_dop853(&n, &fcn, &x, cy, &xend, crtol, catol, &itol, &solout,
        &iout, work, &lwork, iwork, &liwork, crpar, cipar, &idid);
    (*env)->SetDoubleArrayRegion(env, y, 0 , n, cy);
    jclass thisClass = (*env)->GetObjectClass(env, thisObj);
    jfieldID fid = (*env)->GetFieldID(env, thisClass, "y", "[D");
    if (fid == NULL) {
        printf("Field not found.\n");
        return;
    }
    (*env)->SetObjectField(env, thisObj, fid, y);
    return;
}
