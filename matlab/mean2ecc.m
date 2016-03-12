function E = mean2ecc(M, ecc)
keplereq = @(x) x - ecc*sin(x) - M;
keplerderiv = @(x) 1 - ecc*cos(x);
E = newton(M, keplereq, keplerderiv);
