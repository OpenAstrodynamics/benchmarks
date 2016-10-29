package com.helgeeichhorn.icatt;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

public class ElementsFast {
    public static double[] fromRv(Vector3D r, Vector3D v, double mu) {

        final double r_mag = r.getNorm();
        final double v_mag2 = v.getNormSq();
        final Vector3D h = r.crossProduct(v);
        final double h_mag2 = h.getNormSq();
        final double h_mag = FastMath.sqrt(h_mag2);
        final Vector3D k = Vector3D.PLUS_K;
        final Vector3D n = k.crossProduct(h);
        final double n_mag = n.getNorm();
        final double xi = v_mag2/2d - mu/r_mag;
        final Vector3D e = new Vector3D(v_mag2 / mu - 1/r_mag, r, -r.dotProduct(v) / mu, v);
        final double[] ele = new double[6];
        ele[1] = e.getNorm();
        if (ele[1] != 1d) {
            ele[0] = -mu/(2*xi);
        } else {
            ele[0] = h_mag2/mu;
        }
        ele[2] = FastMath.acos(h.getZ()/h_mag);
        ele[3] = FastMath.acos(n.getX()/n_mag);
        ele[4] = FastMath.acos(n.dotProduct(e)/(ele[1]*n_mag));
        ele[5] = FastMath.acos(e.dotProduct(r)/(ele[1]*r_mag));
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
