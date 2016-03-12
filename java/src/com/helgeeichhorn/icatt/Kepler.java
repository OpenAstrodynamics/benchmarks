package com.helgeeichhorn.icatt;

public class Kepler {
    public static double meanMotion(double period, double deltaT) {
        return 2*Math.PI/period*deltaT;
    }

    public static double period(double sma, double mu) {
        return 2*Math.PI*Math.sqrt(Math.pow(sma,3)/mu);
    }

    static class KeplerEquation implements Newton.NewtonInterface {
        private double ecc;
        private double M;
        public KeplerEquation(double meanAnomaly, double eccentricity) {
            M = meanAnomaly;
            ecc = eccentricity;
        }
        public double func(double E) {
            return E - this.ecc*Math.sin(E) - this.M;
        }
        public double deriv(double E) {
            return 1 - this.ecc*Math.cos(E);
        }
    }

    public static double meanToEcc(double M, double ecc) {
        KeplerEquation keplerEq = new KeplerEquation(M, ecc);
        return Newton.getRoot(M, keplerEq);
    }
}
