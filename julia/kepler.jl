export meantoecc, benchmark_kepler

function benchmark_kepler(times::Int)
    r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03]
    v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01]
    mu = 3.986004418e5
    el = elements(r, v, mu)
    worst = -Inf
    best = Inf
    total = 0.0
    meantoecc(pi/2, el[2])
    for i = 1:times
        gc_enable(false)
        t = @elapsed meantoecc(pi/2, el[2])
        gc_enable(true)
        if t > worst
            worst = t
        end
        if t < best
            best = t
        end
        total += t
    end
    println("[$(total/times),$best,$worst]")
end

function period(a::Float64, mu::Float64)
    return 2pi*sqrt(a^3/mu)
end

function newton(x0::Float64, func::Function, derivative::Function, maxiter::Int=50, tol::Float64=sqrt(eps()))
    p0 = x0
    for i = 1:maxiter
        p = p0 - func(p0)/derivative(p0)
        if abs(p - p0) < tol
            return p
        end
        p0 = p
    end
    error("Not converged.")
end

function meantoecc(M::Float64, ecc::Float64)
    kepler(E::Float64) = E - ecc*sin(E) - M
    kepler_der(E::Float64) = 1 - ecc*cos(E)
    return newton(M, kepler, kepler_der)
end
