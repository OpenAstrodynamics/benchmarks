module kepler

use newton
use elements, only: fromrv

implicit none

double precision, parameter :: pi = 3.1415926535897931d0

type, extends(newtoncallback) :: keplercb
    double precision :: ecc
    double precision :: mean
contains
    procedure :: func => keplereq
    procedure :: deriv => keplerderiv
end type keplercb

contains

double precision function mean2ecc(mean, ecc)
    double precision, intent(in) :: mean
    double precision, intent(in) :: ecc

    type(keplercb) :: cb

    cb = keplercb(ecc, mean)
    mean2ecc = newtonsolver(cb, mean)
end function mean2ecc

double precision function keplereq(this, x)
    class(keplercb), intent(in) :: this
    double precision, intent(in) :: x
    keplereq = x - this%ecc * sin(x) - this%mean
end function keplereq

double precision function keplerderiv(this, x)
    class(keplercb), intent(in) :: this
    double precision, intent(in) :: x
    keplerderiv = 1 - this%ecc * cos(x)
end function keplerderiv

double precision function period(sma, mu)
    double precision, intent(in) :: sma
    double precision, intent(in) :: mu
    period = 2*pi*sqrt(sma**3/mu);
end function period

subroutine benchmark_kepler(times)
    integer, intent(in) :: times

    double precision, dimension(3) :: r
    double precision, dimension(3) :: v
    double precision, dimension(6) :: el
    double precision, parameter :: mu = 3.986004418d5
    double precision :: e

    double precision :: start, finish, current
    double precision :: total, best, worst
    integer :: i

    r = [8.59072560d+02, -4.13720368d+03, 5.29556871d+03]
    v = [7.37289205d+00, 2.08223573d+00, 4.39999794d-01]
    el = fromrv(r, v, mu)

    worst = -1d20
    best = 1d20
    total = 0d0
    do i=1, times
        call cpu_time(start)

        e = mean2ecc(pi/2, el(2))

        call cpu_time(finish)
        current = finish - start
        if (current < best .and. current > 0) then
            best = current
        end if
        if (current > worst) then
            worst = current
        end if
        total = total + current
    end do
    print *, total/times, best, worst
end subroutine benchmark_kepler

end module kepler
