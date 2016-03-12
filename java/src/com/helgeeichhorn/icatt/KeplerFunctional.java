package com.helgeeichhorn.icatt;

import java.util.function.Function;

public class KeplerFunctional {
    public static double meanToEcc(double M, double ecc) {
        Function<Double,Double> keplerEq = E -> E - ecc * Math.sin(E) - M;
        Function<Double,Double> keplerDeriv = E -> 1 - ecc * Math.cos(E);
        return NewtonFunctional.getRoot(M, keplerEq, keplerDeriv);
    }
}
