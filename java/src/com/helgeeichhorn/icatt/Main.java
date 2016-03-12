package com.helgeeichhorn.icatt;

public class Main {
    public static void main(String[] args) {
        int n = 100000;
        Benchmark b = new Benchmark();
        b.benchmarkElements(n);
        b.benchmarkKepler(n);
        b.benchmarkKeplerFunctional(n);
        b.benchmarkDopri(n);
        b.benchmarkLambert(n);
    }
}
