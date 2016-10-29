function f = gravity(y, mu)
    r = norm(y(1:3));
    f(1:3) = y(4:6);
    f(4:6) = -mu * y(1:3) / r^3;
end
