using LinearAlgebra
using SpecialFunctions: gamma

export lambert, benchmark_lambert

function benchmark_lambert(times::Int)
    r0 = [5000.0, 10000.0, 2100.0]
    r = [-14600.0, 2500.0, 7000.0]
    tof = 3600.0
    mu = 3.986004418e5
    worst = -Inf
    best = Inf
    total = 0.0
    lambert(mu, r0, r, tof)
    for i = 1:times
        GC.enable(false)
        t = @elapsed lambert(mu, r0, r, tof)
        GC.enable(true)
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

function lambert(k::Float64, r0::Vector{Float64}, r::Vector{Float64}, tof::Float64, short::Bool=true, numiter::Int=35, rtol::Float64=1e-8)
    if short
        t_m = 1
    else
        t_m = -1
    end

    norm_r0 = norm(r0)
    norm_r = norm(r)
    cos_dnu = dot(r0, r) / (norm_r0 * norm_r)

    A = t_m * sqrt(norm_r * norm_r0 * (1 + cos_dnu))

    if A == 0.0
        error("Cannot compute orbit, phase angle is 180 degrees")
    end

    psi = 0.0
    psi_low = -4*pi
    psi_up = 4*pi

    count = 0
    converged = false
    y = 0.0
    while count < numiter
        y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi))
        if A > 0.0 && y < 0.0
            while y < 0.0
                psi_low = psi
                psi = (0.8 * (1.0 / c3(psi)) *
                    (1.0 - (norm_r0 + norm_r) * sqrt(c2(psi)) / A))
                y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi))
            end
        end
        xi = sqrt(y / c2(psi))
        tof_new = (xi^3 * c3(psi) + A * sqrt(y)) / sqrt(k)

        if abs((tof_new - tof) / tof) < rtol
            converged = true
            break
        else
            count += 1
            if tof_new <= tof
                psi_low = psi
            else
                psi_up = psi
            end
            psi = (psi_up + psi_low) / 2
        end
    end

    if !converged
        error("Maximum number of iterations reached")
    end

    f = 1 - y / norm_r0
    g = A * sqrt(y / k)
    gdot = 1 - y / norm_r

    v0 = (r - f * r0) / g
    v = (gdot * r - r0) / g
    return v0, v
end

function c2(psi::Float64)
    eps = 1.0
    if psi > eps
        res = (1 - cos(sqrt(psi))) / psi
    elseif psi < -eps
        res = (cosh(sqrt(-psi)) - 1) / (-psi)
    else
        res = 1.0 / 2.0
        delta = (-psi) / gamma(2 + 2 + 1)
        k = 1
        while res + delta != res
            res += delta
            k += 1
            delta = (-psi)^k / gamma(2*k + 2 + 1)
        end
    end
    return res
end

function c3(psi::Float64)
    eps = 1.0
    if psi > eps
        res = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi))
    elseif psi < -eps
        res = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi))
    else
        res = 1.0 / 6.0
        delta = (-psi) / gamma(2 + 3 + 1)
        k = 1
        while res + delta != res
            res += delta
            k += 1
            delta = (-psi)^k / gamma(2*k + 3 + 1)
        end
    end
    return res
end
