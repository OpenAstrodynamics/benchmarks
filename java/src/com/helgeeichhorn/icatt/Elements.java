package com.helgeeichhorn.icatt;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

public class Elements {
    public static double[] fromRv(Vector3D r, Vector3D v, double mu) {
        Vector3D h, k, n, e;
        double r_mag, v_mag2, h_mag, n_mag, xi;
        double[] ele;

        r_mag = r.getNorm();
        v_mag2 = Math.pow(v.getNorm(), 2d);
        h = r.crossProduct(v);
        h_mag = h.getNorm();
        k = new Vector3D(0,0,1d);
        n = k.crossProduct(h);
        n_mag = n.getNorm();
        xi = v_mag2/2d - mu/r_mag;
        e = r.scalarMultiply(v_mag2 - mu/r_mag).subtract(r.dotProduct(v), v).scalarMultiply(1/mu);
        ele = new double[6];
        ele[1] = e.getNorm();
        if (ele[1] != 1d) {
            ele[0] = -mu/(2*xi);
        } else {
            ele[0] = Math.pow(h_mag,2)/mu;
        }
        ele[2] = Math.acos(h.getZ()/h_mag);
        ele[3] = Math.acos(n.getX()/n_mag);
        ele[4] = Math.acos(n.dotProduct(e)/(ele[1]*n_mag));
        ele[5] = Math.acos(e.dotProduct(r)/(ele[1]*r_mag));
        if (n.getY() < 0) {
            ele[3] = 2*Math.PI - ele[3];
        }
        if (e.getZ() < 0) {
            ele[4] = 2*Math.PI - ele[4];
        }
        if (r.dotProduct(v) < 0) {
            ele[5] = 2*Math.PI - ele[5];
        }
        return ele;
    }
}
