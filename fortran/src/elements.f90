module elements

implicit none

double precision, parameter :: pi = 3.1415926535897931d0

contains

pure function fromrv(r, v, mu) result(ele)
    double precision, dimension(:), intent(in) :: r
    double precision, dimension(:), intent(in) :: v
    double precision, intent(in) :: mu
    double precision, dimension(6) :: ele

    double precision :: r_mag, v_mag, h_mag, n_mag, xi
    double precision, dimension(3) :: h, n, k, e

    r_mag = norm2(r)
    v_mag = norm2(v)
    h = cross(r,v)
    h_mag = norm2(h)
    k = [0d0, 0d0, 1d0]
    n = cross(k, h)
    n_mag = norm2(n)
    xi = v_mag**2/2 - mu/r_mag
    e = ((v_mag**2 - mu/r_mag)*r - v*dot_product(r,v))/mu
    ele(2) = norm2(e)
    if (ele(2) /= 1.0) then
        ele(1) = -mu/(2*xi)
    else
        ele(1) = h_mag**2/mu
    end if
    ele(3) = acos(h(3)/h_mag)
    ele(4) = acos(n(1)/n_mag)
    ele(5) = acos(dot_product(n,e)/(ele(2)*n_mag))
    ele(6) = acos(dot_product(e,r)/(ele(2)*r_mag))
    if (n(2) < 0) then
        ele(4) = 2*pi - ele(4)
    end if
    if (e(3) < 0) then
        ele(5) = 2*pi - ele(5)
    end if
    if (dot_product(r,v) < 0) then
        ele(6) = 2*pi - ele(6)
    end if
end function fromrv

pure function cross(a, b)
    double precision, dimension(:), intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, dimension(3) :: cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

subroutine benchmark_elements(times)
    integer, intent(in) :: times

    double precision, dimension(3) :: r
    double precision, dimension(3) :: v
    double precision, dimension(6) :: el
    double precision, parameter :: mu = 3.986004418d5

    double precision :: current, rate
    integer(kind=8) :: start, finish, irate
    double precision :: total, best, worst
    integer :: i

    r = [8.59072560d+02, -4.13720368d+03, 5.29556871d+03]
    v = [7.37289205d+00, 2.08223573d+00, 4.39999794d-01]

    worst = -1d20
    best = 1d20
    total = 0d0
    call system_clock(count_rate=irate)
    rate = dble(irate)
    do i=1, times
        call system_clock(start)

        el = fromrv(r, v, mu)

        call system_clock(finish)
        current = (finish - start)/rate
        if (current < best .and. current > 0) then
            best = current
        end if
        if (current > worst) then
            worst = current
        end if
        total = total + current
    end do
    print *, total/times, best, worst
end subroutine benchmark_elements

end module elements
