# Account for non-standard library and code locations
push!(LOAD_PATH, abspath(splitdir(@__FILE__)[1]))
push!(Libdl.DL_LOAD_PATH, abspath(joinpath(splitdir(@__FILE__)[1], "..", "lib", "dopri", "build")))

using ICATT

n = 100000
benchmark_elements(n)
benchmark_kepler(n)
benchmark_lambert(n)
benchmark_dopri(n)
