function res = c2(psi)
eps = 1.0;
if psi > eps
    res = (1 - cos(sqrt(psi))) / psi;
elseif psi < -eps
    res = (cosh(sqrt(-psi)) - 1) / (-psi);
else
    res = 1.0 / 2.0;
    delta = (-psi) / gamma(2 + 2 + 1);
    k = 1;
    while res + delta ~= res
        res = res + delta;
        k = k + 1;
        delta = (-psi)^k / gamma(2*k + 2 + 1);
    end
end