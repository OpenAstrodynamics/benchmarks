export dopri, benchmark_dopri

function benchmark_dopri(times::Int)
    y = [8.59072560e+02, -4.13720368e+03, 5.29556871e+03, 7.37289205e+00, 2.08223573e+00, 4.39999794e-01]
    y0 = copy(y)
    mu = 3.986004418e5
    el = elements(y[1:3], y[4:6], mu)
    rpar = collect(mu)
    x = 0.0
    xend = period(el[1], mu)
    worst = -Inf
    best = Inf
    total = 0.0
    dopri!(gravity!, x, y, xend, rpar)
    for i = 1:times
        gc_enable(false)
        t = @elapsed dopri!(gravity!, x, y, xend, rpar)
        gc_enable(true)
        y = copy(y0)
        x = 0.0
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

function gravity!(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint})
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    rpar = unsafe_load(_rpar, 1)
    y = pointer_to_array(_y, n)
    f = pointer_to_array(_f, n)

    r = sqrt(y[1]*y[1]+y[2]*y[2]+y[3]*y[3])
    r3 = r*r*r
    f[1] = y[4]
    f[2] = y[5]
    f[3] = y[6]
    f[4] = -rpar[1] * y[1] / r3
    f[5] = -rpar[1] * y[2] / r3
    f[6] = -rpar[1] * y[3] / r3
    return nothing
end

function _solout(_nr::Ptr{Cint}, _xold::Ptr{Cdouble}, _x::Ptr{Cdouble},
    _y::Ptr{Cdouble}, _n::Ptr{Cint}, _con::Ptr{Cdouble}, _icomp::Ptr{Cint},
    _nd::Ptr{Cint}, _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint}, _irtrn::Ptr{Cint},
    _xout::Ptr{Cdouble})
    return nothing
end

cfcn = cfunction(gravity!, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}))
csolout = cfunction(_solout, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
    Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint},
    Ptr{Cdouble}))

function dopri!(func::Function, x::Float64, y::Vector{Float64}, xend::Float64,
    rpar::Vector{Float64}=Float64[], reltol::Float64=1e-6, abstol::Float64=1e-8)
    n = length(y)
    lwork = 11*n + 8*n + 21
    liwork = n + 21
    work = zeros(Cdouble, lwork)
    iwork = zeros(Cint, liwork)
    ipar = Cint[]
    iout = 0
    idid = 0
    itol = 0
    rtol = collect(reltol)
    atol = collect(abstol)
    ccall((:c_dop853, :libdopri), Void, (Ptr{Cint}, Ptr{Void}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Void}, Ptr{Cint}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Cint}),
    &n, cfcn, &x, y, &xend, rtol, atol, &itol, csolout, &iout,
    work, &lwork, iwork, &liwork, rpar, ipar, &idid)
end
