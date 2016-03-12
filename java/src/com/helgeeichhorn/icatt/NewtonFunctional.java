package com.helgeeichhorn.icatt;

import java.util.function.Function;

public class NewtonFunctional {
    public static double getRoot(
            double p0, Function<Double,Double> func, Function<Double,Double> deriv, int maxiter, double tol) {
        Double result = Double.NaN;
        for (int i=0; i < maxiter; i++) {
            double p = p0 - func.apply(p0) / deriv.apply(p0);
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

    public static double getRoot(double x0, Function<Double,Double> func, Function<Double,Double> deriv) {
        return getRoot(x0, func, deriv, 50, 1e-8);
    }

}
