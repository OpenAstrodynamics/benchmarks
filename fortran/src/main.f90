program main

use elements, only: benchmark_elements
use kepler, only: benchmark_kepler
use fdopri, only: benchmark_dopri
use lambert, only: benchmark_lambert

implicit none

integer, parameter :: times = 100000

call benchmark_elements(times)
call benchmark_kepler(times)
call benchmark_lambert(times)
call benchmark_dopri(times)

end program main
