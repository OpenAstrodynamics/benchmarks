module dopri

use iso_c_binding, only: c_double, c_int, c_funptr, c_f_procpointer, c_associated, c_loc

implicit none

abstract interface
    subroutine dopfcn(n, x, y, f, rpar, ipar)
        integer, intent(inout) :: n
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
        double precision, dimension(:), intent(inout) :: rpar
        double precision, dimension(n), intent(inout) :: f
        integer, dimension(:), intent(inout) :: ipar
    end subroutine dopfcn
    subroutine dopsolout(nr, xold, x, y, n, con, icomp,&
            nd, rpar, ipar, irtrn, xout)
        integer, intent(inout) :: nr
        double precision, intent(inout) :: xold
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
        integer, intent(inout) :: n
        double precision, dimension(8*nd), intent(inout) :: con
        integer, dimension(nd), intent(inout) :: icomp
        integer, intent(inout) :: nd
        double precision, dimension(:), intent(inout) :: rpar
        integer, dimension(:), intent(inout) :: ipar
        integer, intent(inout) :: irtrn
        double precision, intent(inout) :: xout
    end subroutine dopsolout
end interface

interface
    subroutine dop853(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
        procedure(dopfcn), pointer :: fcn
        procedure(dopsolout), pointer :: solout
        integer, intent(inout) :: n
        integer, intent(inout) :: itol
        integer, intent(inout) :: iout
        integer, intent(inout) :: lwork
        integer, intent(inout) :: liwork
        integer, dimension(:),intent(inout) :: ipar
        integer, intent(inout) :: idid
        double precision, intent(inout) :: xend
        double precision, dimension(n), intent(inout) :: rtol
        double precision, dimension(n), intent(inout) :: atol
        double precision, dimension(:), intent(inout) :: rpar
        double precision, dimension(lwork), intent(inout) :: work
        integer, dimension(liwork), intent(inout) :: iwork
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
    end subroutine dop853

    double precision function contd8(ii, x, con, icomp, nd)
        integer, intent(inout) :: ii
        double precision, intent(inout) :: x
        double precision, dimension(8*nd), intent(inout) :: con
        integer, dimension(nd), intent(inout) :: icomp
        integer, intent(inout) :: nd
    end function contd8
end interface

abstract interface
    subroutine c_fcn(n, x, y, f, rpar, ipar)
        import :: c_int
        import :: c_double
        integer(c_int), intent(inout) :: n
        integer(c_int), dimension(:),intent(inout) :: ipar
        real(c_double), intent(inout) :: x
        real(c_double), dimension(n), intent(inout) :: y
        real(c_double), dimension(:), intent(inout) :: rpar
        real(c_double), dimension(n), intent(inout) :: f
    end subroutine c_fcn
    subroutine c_solout(nr, xold, x, y, n, con, icomp,&
            nd, rpar, ipar, irtrn, xout)
        import :: c_int
        import :: c_double
        integer(c_int), intent(inout) :: n
        integer(c_int), intent(inout) :: nr
        integer(c_int), intent(inout) :: nd
        integer(c_int), intent(inout) :: irtrn
        integer(c_int), dimension(:), intent(inout) :: ipar
        integer(c_int), dimension(nd), intent(inout) :: icomp
        real(c_double), intent(inout) :: xold
        real(c_double), intent(inout) :: x
        real(c_double), dimension(n), intent(inout) :: y
        real(c_double), dimension(8*nd), intent(inout) :: con
        real(c_double), dimension(:), intent(inout) :: rpar
        real(c_double), intent(inout) :: xout
    end subroutine c_solout
end interface

contains

subroutine c_dop853(n, cfcn, x, y, xend, rtol, atol,&
        itol, csolout, iout, work, lwork, iwork,&
        liwork, rpar, ipar, idid) bind(c)
    integer(c_int), intent(inout) :: n
    type(c_funptr), intent(in), value :: cfcn
    real(c_double), intent(inout) :: x
    real(c_double), dimension(n), intent(inout) :: y
    real(c_double), intent(inout) :: xend
    real(c_double), dimension(n), intent(inout) :: rtol
    real(c_double), dimension(n), intent(inout) :: atol
    integer(c_int), intent(inout) :: itol
    type(c_funptr), intent(in), value :: csolout
    integer(c_int), intent(inout) :: iout
    real(c_double), dimension(lwork), intent(inout) :: work
    integer(c_int), intent(inout) :: lwork
    integer(c_int), dimension(liwork), intent(inout) :: iwork
    integer(c_int), intent(inout) :: liwork
    real(c_double), dimension(:), intent(inout) :: rpar
    integer(c_int), dimension(:),intent(inout) :: ipar
    integer(c_int), intent(inout) :: idid

    procedure(c_fcn), pointer :: fcn
    procedure(c_solout), pointer :: solout

    call c_f_procpointer(cfcn, fcn)
    call c_f_procpointer(csolout, solout)

    print *, c_associated(cfcn)
    print *, loc(cfcn)
    print *, associated(fcn)
    print *, loc(fcn)

    call dop853(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
end subroutine c_dop853

function c_contd8(ii, x, con, icomp, nd) result(ret) bind(c)
    real(c_double) :: ret
    integer(c_int), intent(inout) :: ii
    real(c_double), intent(inout) :: x
    real(c_double), dimension(8*nd), intent(inout) :: con
    integer(c_int), dimension(nd), intent(inout) :: icomp
    integer(c_int), intent(inout) :: nd
    ret = contd8(ii, x, con, icomp, nd)
end function c_contd8

end module
