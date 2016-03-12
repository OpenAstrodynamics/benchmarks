module newton

implicit none

type, abstract :: newtoncallback
contains
    procedure(dpfunc), deferred :: func
    procedure(dpfunc), deferred :: deriv
end type newtoncallback

abstract interface
    function dpfunc(this, x) result(fx)
        import :: newtoncallback
        class(newtoncallback), intent(in) :: this
        double precision, intent(in) :: x
        double precision :: fx
    end function dpfunc
end interface

contains

double precision function newtonsolver(cb, x0, maxiter, tol)
    class(newtoncallback), intent(in) :: cb
    double precision, intent(in) :: x0
    integer, intent(in), optional :: maxiter
    double precision, intent(in), optional :: tol

    integer :: maxiter_
    double precision :: tol_

    integer :: i
    double precision :: p
    double precision :: p0

    p0 = x0
    maxiter_ = 50
    if (present(maxiter)) maxiter_ = maxiter
    tol_ = 1d-8
    if (present(tol)) tol_ = tol

    do i = 1, maxiter_
        p = p0 - cb.func(p0) / cb.deriv(p0)
        if (abs(p - p0) < tol_) then
            newtonsolver = p
            return
        end if
        p0 = p
    end do

    write(*,*) "Not converged."
    stop 1
end function newtonsolver

end module newton
