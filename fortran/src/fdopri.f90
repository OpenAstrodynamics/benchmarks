module fdopri

use dopri
use elements, only: fromrv
use kepler, only: period

implicit none

contains

subroutine integrate(fcn, y, t, tend, rpar, ipar,&
        stat, solout, rtol, atol)
    double precision, dimension(:), intent(inout) :: y
    double precision, intent(inout) :: t
    double precision, intent(inout) :: tend
    double precision, dimension(:), intent(inout) :: rpar
    integer, dimension(:), intent(inout) :: ipar
    integer, intent(out), optional :: stat
    double precision, dimension(:), intent(in), optional :: rtol
    double precision, dimension(:), intent(in), optional :: atol
    procedure(dopfcn) :: fcn
    procedure(dopsolout), optional :: solout

    integer :: lwork
    integer :: liwork
    integer :: iout
    integer :: idid
    integer :: itol
    integer, dimension(:), allocatable :: iwork
    integer :: n
    double precision, dimension(:), allocatable :: rtol_
    double precision, dimension(:), allocatable :: atol_
    double precision, dimension(:), allocatable :: work
    double precision, dimension(:), allocatable :: y0

    n = size(y)
    lwork = 11*n + 8*n + 21
    liwork = n + 21
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(y0(n))
    allocate(rtol_(n))
    allocate(atol_(n))

    if (present(rtol)) then
        rtol_ = rtol
    else
        rtol_ = 1d-6
    end if
    if (present(atol)) then
        atol_ = atol
    else
        atol_ = 1d-8
    end if

    iout = 0
    itol = 0
    iwork = 0
    work = 0d0

    y0 = y
    if (present(solout)) then
        call dop853(n, fcn, t, y0, tend, rtol_, atol_,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
    else
        call dop853(n, fcn, t, y0, tend, rtol_, atol_,&
            itol, soldummy, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
    end if
    y = y0

    if (present(stat)) stat = idid
end subroutine integrate

subroutine gravity(n, t, y, f, rpar, ipar)
    integer, intent(inout) :: n
    double precision, intent(inout) :: t
    double precision, dimension(n), intent(inout) :: y
    double precision, dimension(n), intent(inout) :: f
    double precision, dimension(:), intent(inout) :: rpar
    integer, dimension(:),intent(inout) :: ipar

    double precision :: r

    r = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
    f(1) = y(4)
    f(2) = y(5)
    f(3) = y(6)

    f(4) = -rpar(1) * y(1) / r**3
    f(5) = -rpar(1) * y(2) / r**3
    f(6) = -rpar(1) * y(3) / r**3
end subroutine gravity

subroutine soldummy(nr, xold, x, y, n, con, icomp,&
                    nd, rpar, ipar, irtrn, xout)
    integer, intent(inout) :: n
    integer, intent(inout) :: nr
    integer, intent(inout) :: nd
    integer, intent(inout) :: irtrn
    integer, dimension(:), intent(inout) :: ipar
    integer, dimension(nd), intent(inout) :: icomp
    double precision, intent(inout) :: xold
    double precision, intent(inout) :: x
    double precision, dimension(n), intent(inout) :: y
    double precision, dimension(8*nd), intent(inout) :: con
    double precision, dimension(:), intent(inout) :: rpar
    double precision, intent(inout) :: xout
    xout = 0d0
end subroutine soldummy

subroutine benchmark_dopri(times)
    integer, intent(in) :: times

    double precision, dimension(6) :: rv
    double precision, dimension(6) :: rv0
    double precision, dimension(6) :: el
    double precision, dimension(1) :: rpar
    integer, dimension(1) :: ipar
    double precision :: t
    double precision :: tend
    double precision, parameter :: mu = 3.986004418d5

    double precision :: current, rate
    integer(kind=8) :: start, finish, irate
    double precision :: total, best, worst
    integer :: i

    rv = [8.59072560d+02, -4.13720368d+03, 5.29556871d+03, 7.37289205d+00, 2.08223573d+00, 4.39999794d-01]
    rv0 = rv
    el = fromrv(rv(:3), rv(4:), mu)
    rpar(1) = mu
    ipar(1) = 0
    t = 0d0
    tend = period(el(1), mu)

    worst = -1d20
    best = 1d20
    total = 0d0
    call system_clock(count_rate=irate)
    rate = dble(irate)
    do i=1, times
        call system_clock(start)

        call integrate(gravity, rv, t, tend, rpar, ipar)

        call system_clock(finish)
        current = (finish - start)/rate
        rv = rv0
        t = 0d0
        if (current < best .and. current > 0) then
            best = current
        end if
        if (current > worst) then
            worst = current
        end if
        total = total + current
    end do
    print *, total/times, best, worst
end subroutine benchmark_dopri

end module fdopri
