package com.helgeeichhorn.icatt;

public class Newton {
    public interface NewtonInterface {
        double func(double x);
        double deriv(double x);
    }

    public static double getRoot(double p0, NewtonInterface newton, int maxiter, double tol) {
        Double result = Double.NaN;
        for (int i=0; i < maxiter; i++) {
            double p = p0 - newton.func(p0) / newton.deriv(p0);
            if (Math.abs(p - p0) < tol) {
                result = p;
                break;
            }
            p0 = p;
        }
        if (result.isNaN()) {
            throw new RuntimeException("Not converged.");
        } else {
            return result;
        }
    }

    public static double getRoot(double x0, NewtonInterface newton) {
       return getRoot(x0, newton, 50, 1e-8);
    }
}
