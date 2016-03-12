module lambert

implicit none

double precision, parameter :: pi = 3.1415926535897931d0

type lambertresult
    double precision, dimension(3) :: v0
    double precision, dimension(3) :: v
end type lambertresult

contains

type(lambertresult) function lambertsolver(k, r0, r, tof, short, numiter, rtol)
    double precision, intent(in) :: k
    double precision, dimension(:), intent(in) :: r0
    double precision, dimension(:), intent(in) :: r
    double precision, intent(in) :: tof
    logical, intent(in), optional :: short
    integer, intent(in), optional :: numiter
    double precision, intent(in), optional :: rtol

    logical :: short_
    integer :: numiter_
    double precision :: rtol_

    double precision :: norm_r0, norm_r, cos_dnu
    double precision :: a, psi, psi_low, psi_up
    double precision :: y, xi, tof_new, g, gdot
    integer :: t_m, counter
    double precision, dimension(3) :: f, v, v0

    short_ = .true.
    if (present(short)) short_ = short
    numiter_ = 35
    if (present(numiter)) numiter_ = numiter
    rtol_ = 1e-8
    if (present(rtol)) rtol_ = rtol

    if (short_) then
        t_m = 1
    else
        t_m = -1
    end if

    norm_r0 = norm2(r0)
    norm_r = norm2(r)
    cos_dnu = dot_product(r0, r) / (norm_r0 * norm_r)

    a = t_m * sqrt(norm_r * norm_r0 * (1 + cos_dnu))

    if (a == 0d0) then
        write(*,*) "Cannot compute orbit, phase angle is 180 degrees"
        stop 1
    end if

    psi = 0d0
    psi_low = -4 * pi
    psi_up = 4 * pi

    counter = 0
    do while (counter < numiter_)
        y = norm_r0 + norm_r + a * (psi * c3(psi) - 1) / sqrt(c2(psi))
        if (a > 0d0 .and. y < 0d0) then
            do while (y < 0d0)
                psi_low = psi
                psi = (0.8d0 * (1d0 / c3(psi)) * &
                    (1d0 - (norm_r0 + norm_r) * sqrt(c2(psi)) / a))
                y = norm_r0 + norm_r + a * (psi * c3(psi) - 1) / sqrt(c2(psi))
            end do
        end if
        xi = sqrt(y / c2(psi))
        tof_new = (xi**3 * c3(psi) + A * sqrt(y)) / sqrt(k)

        if (abs((tof_new - tof) / tof) < rtol_) then
            exit
        else
            counter = counter + 1
            if (tof_new <= tof) then
                psi_low = psi
            else
                psi_up = psi
            end if
            psi = (psi_up + psi_low) / 2
        end if
    end do

    if (counter > numiter_) then
        write(*,*) "Maximum number of iterations reached."
        stop 1
    end if

     f = 1 - y / norm_r0
     g = a * sqrt(y / k)
     gdot = 1 - y / norm_r

     v0 = (r - f * r0) / g
     v = (gdot * r - r0) / g
     lambertsolver = lambertresult(v0, v)
end function lambertsolver

double precision function c2(psi)
    double precision, intent(in) :: psi

    double precision :: eps
    double precision :: delta
    integer :: k

    eps = 1.0
    if (psi > eps) then
        c2 = (1 - cos(sqrt(psi))) / psi
    else if (psi < -eps) then
        c2 = (cosh(sqrt(-psi)) - 1) / (-psi)
    else
        c2 = 1d0 / 2d0
        delta = (-psi) / gamma(2d0 + 2d0 + 1d0)
        k = 1
        do while (c2 + delta /= c2)
            c2 = c2 + delta
            k = k + 1
            delta = (-psi)**k / gamma(2 * k + 2d0 + 1d0)
        end do
    end if
end function c2

double precision function c3(psi)
    double precision, intent(in) :: psi

    double precision :: eps
    double precision :: delta
    integer :: k

    eps = 1.0
    if (psi > eps) then
        c3 = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi))
    else if (psi < -eps) then
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi))
    else
        c3 = 1d0 / 6d0
        delta = (-psi) / gamma(2d0 + 3d0 + 1d0)
        k = 1
        do while (c3 + delta /= c3)
            c3 = c3 + delta
            k = k + 1
            delta = (-psi) ** k / gamma(2 * k + 3d0 + 1d0)
        end do
    end if
end function c3

subroutine benchmark_lambert(times)
    integer, intent(in) :: times

    double precision, dimension(3) :: r
    double precision, dimension(3) :: r0
    double precision :: tof
    double precision, parameter :: mu = 3.986004418d5
    type(lambertresult) :: res

    double precision :: start, finish, current
    double precision :: total, best, worst
    integer :: i

    r0 = [5000d0, 10000d0, 2100d0]
    r = [-14600d0, 2500d0, 7000d0]
    tof = 3600d0

    worst = -1d20
    best = 1d20
    total = 0d0
    do i=1, times
        call cpu_time(start)

        res = lambertsolver(mu, r0, r, tof)

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
end subroutine benchmark_lambert

end module lambert
