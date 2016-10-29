function benchmark(times)

r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03];
v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01];
mu = 3.986004418e5;
rlam0 = [5000.0, 10000.0, 2100.0];
rlam = [-14600.0, 2500.0, 7000.0];
tof = 3600.0;

el = elements(r, v, mu);
tend = period(el(1), mu)

worst = -inf;
best = inf;
total = 0;
for ii = 1:times
    tic;
    elements(r, v, mu);
    t = toc;

    total = total + t;
    if t > worst
        worst = t;
    end
    if t < best
        best = t;
    end
end

disp(['[',num2str(total/times),',',num2str(best),',',num2str(worst),']'])

worst = -inf;
best = inf;
total = 0;
for ii = 1:times
    tic;
    mean2ecc(pi/2, el(2));
    t = toc;

    total = total + t;
    if t > worst
        worst = t;
    end
    if t < best
        best = t;
    end
end

disp(['[',num2str(total/times),',',num2str(best),',',num2str(worst),']'])

worst = -inf;
best = inf;
total = 0;
for ii = 1:times
    tic;
    propagator('gravity', [r,v], 0, tend, mu);
    t = toc;

    total = total + t;
    if t > worst
        worst = t;
    end
    if t < best
        best = t;
    end
end

disp(['[',num2str(total/times),',',num2str(best),',',num2str(worst),']'])

