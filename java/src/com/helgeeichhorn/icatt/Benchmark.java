package com.helgeeichhorn.icatt;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.Arrays;

public class Benchmark {
    private static double mu = 3.986004418e5;
    private static double rv[] = {
            8.59072560e+02, -4.13720368e+03, 5.29556871e+03,
            7.37289205e+00, 2.08223573e+00, 4.39999794e-01};
    public static Vector3D r = new Vector3D(8.59072560e+02, -4.13720368e+03, 5.29556871e+03);
    public static Vector3D v = new Vector3D(7.37289205e+00, 2.08223573e+00, 4.39999794e-01);
    public static void benchmarkElements(int times) {
        double worst = Double.NEGATIVE_INFINITY;
        double best = Double.POSITIVE_INFINITY;
        double all = 0;
        Long current;
        double seconds;
        long start;
        long end;
        ElementsFast.fromRv(r, v, mu);
        for (int i = 0; i < times; i++) {
            start = System.nanoTime();
            ElementsFast.fromRv(r, v, mu);
            end = System.nanoTime();
            current = end - start;
            seconds = current.doubleValue()/1e9;
            if (seconds < best) {
                best = seconds;
            }
            if (seconds > worst) {
                worst = seconds;
            }
            all += seconds;

        }
        double[] results = {all/times, best, worst};
        System.out.println(Arrays.toString(results));
    }
    public static void benchmarkKepler(int times) {
        double[] ele = Elements.fromRv(r, v, mu);
        double worst = Double.NEGATIVE_INFINITY;
        double best = Double.POSITIVE_INFINITY;
        double all = 0;
        Long current;
        double seconds;
        long start;
        long end ;
        Kepler.meanToEcc(Math.PI/2, ele[1]);
        for (int i = 0; i < times; i++) {
            start = System.nanoTime();
            Kepler.meanToEcc(Math.PI/2, ele[1]);
            end = System.nanoTime();
            current = end - start;
            seconds = current.doubleValue()/1e9;
            if (seconds < best) {
                best = seconds;
            }
            if (seconds > worst) {
                worst = seconds;
            }
            all += seconds;

        }
        double[] results = {all/times, best, worst};
        System.out.println(Arrays.toString(results));
    }
    public static void benchmarkKeplerFunctional(int times) {
        double[] ele = Elements.fromRv(r, v, mu);
        double worst = Double.NEGATIVE_INFINITY;
        double best = Double.POSITIVE_INFINITY;
        double all = 0;
        Long current;
        double seconds;
        long start;
        long end ;
        KeplerFunctional.meanToEcc(Math.PI/2, ele[1]);
        for (int i = 0; i < times; i++) {
            start = System.nanoTime();
            KeplerFunctional.meanToEcc(Math.PI/2, ele[1]);
            end = System.nanoTime();
            current = end - start;
            seconds = current.doubleValue()/1e9;
            if (seconds < best) {
                best = seconds;
            }
            if (seconds > worst) {
                worst = seconds;
            }
            all += seconds;

        }
        double[] results = {all/times, best, worst};
        System.out.println(Arrays.toString(results));
    }
    public static void benchmarkDopri(int times) {
        double[] ele = Elements.fromRv(r, v, mu);
        double period = Kepler.period(ele[0], mu);
        double worst = Double.NEGATIVE_INFINITY;
        double best = Double.POSITIVE_INFINITY;
        double all = 0;
        Long current;
        double seconds;
        long start;
        long end ;
        double rpar[] = {mu};
        Dopri d = new Dopri(new Gravity(), 0, rv, rpar);
        d.integrate(period);
        for (int i = 0; i < times; i++) {
            start = System.nanoTime();
            d = new Dopri(new Gravity(), 0, rv, rpar);
            d.integrate(period);
            end = System.nanoTime();
            current = end - start;
            seconds = current.doubleValue()/1e9;
            if (seconds < best) {
                best = seconds;
            }
            if (seconds > worst) {
                worst = seconds;
            }
            all += seconds;

        }
        double[] results = {all/times, best, worst};
        System.out.println(Arrays.toString(results));
    }
    public static void benchmarkLambert(int times) {
        Vector3D rlam0 = new Vector3D(5000.0, 10000.0, 2100.0);
        Vector3D rlam = new Vector3D(-14600.0, 2500.0, 7000.0);
        double tof = 3600;
        double worst = Double.NEGATIVE_INFINITY;
        double best = Double.POSITIVE_INFINITY;
        double all = 0;
        Long current;
        double seconds;
        long start;
        long end ;
        LambertFast.solve(mu, rlam0, rlam, tof);
        for (int i = 0; i < times; i++) {
            start = System.nanoTime();
            LambertFast.solve(mu, rlam0, rlam, tof);
            end = System.nanoTime();
            current = end - start;
            seconds = current.doubleValue()/1e9;
            if (seconds < best) {
                best = seconds;
            }
            if (seconds > worst) {
                worst = seconds;
            }
            all += seconds;

        }
        double[] results = {all/times, best, worst};
        System.out.println(Arrays.toString(results));
    }
}
