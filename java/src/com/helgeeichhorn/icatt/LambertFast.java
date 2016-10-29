package com.helgeeichhorn.icatt;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

public class LambertFast {
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
        double A = tm * FastMath.sqrt(normR * normR0 * (1 + cosDnu));

        if (A == 0.0) {
            throw new RuntimeException("Cannot compute orbit, phase angle is 180 degrees");
        }
        double psi = 0.0;
        double psiLow = -4*FastMath.PI;
        double psiUp = 4*FastMath.PI;

        double y = 0;
        int count = 0;
        while (count < numiter) {
            y = normR0 + normR + A * (psi * c3(psi) - 1) / FastMath.sqrt(c2(psi));
            if (A > 0 & y < 0) {
                while (y < 0) {
                    psiLow = psi;
                    psi = (0.8 * (1.0 / c3(psi)) *
                            (1.0 - (normR0 + normR) * FastMath.sqrt(c2(psi)) / A));
                    y = normR + normR0 + A * (psi * c3(psi) - 1) / FastMath.sqrt(c2(psi));
                }
            }
            double sss = y / c2(psi);
            double xi = FastMath.sqrt(sss);
            double tofNew = (sss * xi * c3(psi) + A * FastMath.sqrt(y)) / FastMath.sqrt(k);

            if (FastMath.abs((tofNew - tof) / tof) < rtol) {
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
        double g = A * FastMath.sqrt(y / k);
        double gdot = 1 - y / normR;
        Vector3D v0 = new Vector3D(1/g, r, -f/g, r0);
        Vector3D v = new Vector3D(gdot/g, r, -1/g, r0);
        return new LambertResult(v0, v);
    }

    public static double c2(double psi) {
        double res;
        double eps = 1.0;
        if (psi > eps) {
            res = (1 - FastMath.cos(FastMath.sqrt(psi))) / psi;
        } else if (psi < -eps) {
            res = (FastMath.cosh(FastMath.sqrt(-psi)) - 1) / -psi;
        } else {
            res = 1.0 / 2.0;
            double delta = (-psi) / Gamma.gamma(2+2+1);
            int k = 1;
            double mpsik = -psi;
            while (res + delta != res) {
                res += delta;
                k += 1;
                mpsik *= -psi;
                delta = mpsik / Gamma.gamma(2 * k + 2 + 1);
            }
        }
        return res;
    }

    public static double c3(double psi) {
        double res;
        double eps = 1.0;
        if (psi > eps) {
            res = (FastMath.sqrt(psi) - FastMath.sin(FastMath.sqrt(psi))) / (psi * FastMath.sqrt(psi));
        } else if (psi < -eps) {
            res = (FastMath.sinh(FastMath.sqrt(-psi)) - FastMath.sqrt(-psi)) / (-psi * FastMath.sqrt(-psi));
        } else {
            res = 1.0/6.0;
            double delta = (-psi) / Gamma.gamma(2 + 3 + 1);
            int k = 1;
            double mpsik = -psi;
            while (res + delta != res) {
                res += delta;
                k += 1;
                mpsik *= -psi;
                delta = mpsik / Gamma.gamma(2 * k + 3 + 1);
            }
        }
        return res;
    }
}
