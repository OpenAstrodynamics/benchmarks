package com.helgeeichhorn.icatt;

import java.util.Arrays;

public class Dopri {
    static {
        System.loadLibrary("jnilib");
    }
    private DopriInterface dopri;
    private double x;
    private double[] y;
    private double[] rtol;
    private double[] atol;
    private double[] rpar;

    Dopri(DopriInterface f, double x0, double[] y0,
          double[] reltol, double[] abstol, double[] params) {
        dopri = f;
        x = x0;
        y = y0;
        rtol = reltol;
        atol = abstol;
        rpar = params;
    }

    Dopri(DopriInterface f, double x0, double[] y0, double[] params) {
        dopri = f;
        x = x0;
        y = y0;
        rtol = new double[1];
        Arrays.fill(rtol, 1e-6);
        atol = new double[1];
        Arrays.fill(atol, 1e-8);
        rpar = params;
    }

    public double[] integrate(double xend) {
        int n = this.y.length;
        int itol = 0;
        int iout = 0;
        int idid = 0;
        int[] ipar = {this.rpar.length};
        jdop853(this.dopri, n, this.x, this.y, xend, this.rtol, this.atol, itol, iout,
                this.rpar, ipar, idid);
        return this.y;
    }

    public native void jdop853(DopriInterface dopri, int n, double x, double[] y, double xend, double[] rtol,
        double[] atol, int itol, int iout, double[] rpar, int[] ipar, int idid);


    public interface DopriInterface {
        double[] func(double x, double[] y, double[] params);
    }
}
