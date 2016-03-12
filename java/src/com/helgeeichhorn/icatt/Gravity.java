package com.helgeeichhorn.icatt;

public class Gravity implements Dopri.DopriInterface {
    public double[] func(double x, double[] y, double[] rpar) {
        double r = Math.sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
        double[] f = new double[6];
        f[0] = y[3];
        f[1] = y[4];
        f[2] = y[5];
        f[3] = -rpar[0]*y[0]/Math.pow(r,3);
        f[4] = -rpar[0]*y[1]/Math.pow(r,3);
        f[5] = -rpar[0]*y[2]/Math.pow(r,3);
        return f;
    }
}
