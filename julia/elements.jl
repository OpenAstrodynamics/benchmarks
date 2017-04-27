using BenchmarkTools
export elements, benchmark_elements

function benchmark_elements(times)
    r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03]
    v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01]
    mu = 3.986004418e5
    elements(r, v, mu)
    b = @benchmark elements($r, $v, $mu)
    println(b)
end
#= function benchmark_elements(times::Int) =#
#=     r = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03] =#
#=     v = [7.37289205e+00, 2.08223573e+00, 4.39999794e-01] =#
#=     mu = 3.986004418e5 =#
#=     worst = -Inf =#
#=     best = Inf =#
#=     total = 0.0 =#
#=     elements(r, v, mu) =#
#=     for i = 1:times =#
#=         gc_enable(false) =#
#=         t = @elapsed elements(r, v, mu) =#
#=         gc_enable(true) =#
#=         if t > worst =#
#=             worst = t =#
#=         end =#
#=         if t < best =#
#=             best = t =#
#=         end =#
#=         total += t =#
#=     end =#
#=     println("[$(total/times),$best,$worst]") =#
#= end =#

function elements(r::AbstractArray, v::AbstractArray, mu::Float64)
    rm = norm(r)
    vm = norm(v)
    h = cross(r, v)
    hm = norm(h)
    k = [0.0, 0.0, 1.0]
    n = cross(k, h)
    nm = norm(n)
    xi = vm^2/2 - mu/rm
    ec = ((vm^2 - mu/rm)*r - v*dot(r, v))/mu
    ecc = norm(ec)
    if ecc != 1
        sma = -mu/(2*xi)
    else
        sma = hm^2/mu
    end
    inc = acos(h[3]/hm)
    node = acos(n[1]/nm)
    peri = acos(dot(n, ec)/(ecc*nm))
    ano = acos(dot(ec, r)/(ecc*rm))
    if n[2] < 0
        node = 2*pi - node
    end
    if ec[3] < 0
        peri = 2*pi - peri
    end
    if dot(r, v) < 0
        ano = 2*pi - ano
    end
    return sma, ecc, inc, node, peri, ano
end
