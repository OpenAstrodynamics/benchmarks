r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03];
v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01];
mu = 3.986004418e5;
el = elements(r, v, mu);
tend = period(el(1), mu);
propagator(1)
rv1 = propagator('gravity', [r,v], 0, tend, mu);
propagator(0)

format long;
disp(rv1)
% exit(0)
