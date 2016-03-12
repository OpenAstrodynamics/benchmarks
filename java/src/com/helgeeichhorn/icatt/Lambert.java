package com.helgeeichhorn.icatt;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.special.Gamma;

public class Lambert {
    public static final class LambertResult {
        public Vector3D v0;
        public Vector3D v;

        LambertResult(Vector3D v0, Vector3D v) {
            this.v0 = v0;
            this.v = v;
        }
    }
    public static LambertResult solve(double k, Vector3D r0, Vector3D r, double tof) {
        boolean shortway = true;
        int numiter = 35;
        double rtol = 1e-8;
        return solve(k, r0, r, tof, shortway, numiter, rtol);
    }
    public static LambertResult solve(double k, Vector3D r0, Vector3D r, double tof,
                                 boolean shortway, int numiter, double rtol) {
        int tm;
        if (shortway) {
            tm = 1;
        } else {
            tm = -1;
        }
        double normR0 = r0.getNorm();
        double normR = r.getNorm();
        double cosDnu = r0.dotProduct(r) / (normR0*normR);
        double A = tm * Math.sqrt(normR * normR0 * (1 + cosDnu));

        if (A == 0.0) {
            throw new RuntimeException("Cannot compute orbit, phase angle is 180 degrees");
        }
        double psi = 0.0;
        double psiLow = -4*Math.PI;
        double psiUp = 4*Math.PI;

        double y = 0;
        int count = 0;
        while (count < numiter) {
            y = normR0 + normR + A * (psi * c3(psi) - 1) / Math.sqrt(c2(psi));
            if (A > 0 & y < 0) {
                while (y < 0) {
                    psiLow = psi;
                    psi = (0.8 * (1.0 / c3(psi)) *
                            (1.0 - (normR0 + normR) * Math.sqrt(c2(psi)) / A));
                    y = normR + normR0 + A * (psi * c3(psi) - 1) / Math.sqrt(c2(psi));
                }
            }
            double xi = Math.sqrt(y / c2(psi));
            double tofNew = (Math.pow(xi, 3) * c3(psi) + A * Math.sqrt(y)) / Math.sqrt(k);

            if (Math.abs((tofNew - tof) / tof) < rtol) {
                break;
            } else {
                count += 1;
                if (tofNew <= tof) {
                    psiLow = psi;
                } else {
                    psiUp = psi;
                }
                psi = (psiUp + psiLow) / 2;
            }
        }
        if (count >= numiter) {
            throw new RuntimeException("Maximum number of iterations reached.");
        }
        double f = 1 - y / normR0;
        double g = A * Math.sqrt(y / k);
        double gdot = 1 - y / normR;
        Vector3D v0 = r.subtract(r0.scalarMultiply(f)).scalarMultiply(1/g);
        Vector3D v = r.scalarMultiply(gdot).subtract(r0).scalarMultiply(1/g);
        return new LambertResult(v0, v);
    }

    public static double c2(double psi) {
        double res;
        double eps = 1.0;
        if (psi > eps) {
            res = (1 - Math.cos(Math.sqrt(psi))) / psi;
        } else if (psi < -eps) {
            res = (Math.cosh(Math.sqrt(-psi)) - 1) / -psi;
        } else {
            res = 1.0 / 2.0;
            double delta = (-psi) / Gamma.gamma(2+2+1);
            int k = 1;
            while (res + delta != res) {
                res += delta;
                k += 1;
                delta = Math.pow(-psi, k) / Gamma.gamma(2 * k + 2 + 1);
            }
        }
        return res;
    }

    public static double c3(double psi) {
        double res;
        double eps = 1.0;
        if (psi > eps) {
            res = (Math.sqrt(psi) - Math.sin(Math.sqrt(psi))) / (psi * Math.sqrt(psi));
        } else if (psi < -eps) {
            res = (Math.sinh(Math.sqrt(-psi)) - Math.sqrt(-psi)) / (-psi * Math.sqrt(-psi));
        } else {
            res = 1.0/6.0;
            double delta = (-psi) / Gamma.gamma(2 + 3 + 1);
            int k = 1;
            while (res + delta != res) {
                res += delta;
                k += 1;
                delta = Math.pow(-psi, k) / Gamma.gamma(2 * k + 3 + 1);
            }
        }
        return res;
    }
}
