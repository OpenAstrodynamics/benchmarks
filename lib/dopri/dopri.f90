module dopri

use iso_c_binding, only: c_double, c_int, c_funptr, c_f_procpointer

implicit none

abstract interface
    subroutine dopfcn(n, x, y, f, rpar, ipar)
        integer, intent(inout) :: n
        integer, dimension(:),intent(inout) :: ipar
        double precision, intent(inout) :: x
        double precision, dimension(n), intent(inout) :: y
        double precision, dimension(:), intent(inout) :: rpar
        double precision, dimension(n), intent(inout) :: f
    end subroutine dopfcn
    subroutine dopsolout(nr, xold, x, y, n, con, icomp,&
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
    end subroutine dopsolout
end interface

interface
    subroutine dop853(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
        interface
            subroutine fcn(n, x, y, f, rpar, ipar)
                integer, intent(inout) :: n
                integer, dimension(:),intent(inout) :: ipar
                double precision, intent(inout) :: x
                double precision, dimension(n), intent(inout) :: y
                double precision, dimension(:), intent(inout) :: rpar
                double precision, dimension(n), intent(inout) :: f
            end subroutine fcn
            subroutine solout(nr, xold, x, y, n, con, icomp,&
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
            end subroutine solout
        end interface
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

    subroutine dopri5(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
        interface
            subroutine fcn(n, x, y, f, rpar, ipar)
                integer, intent(inout) :: n
                integer, dimension(:),intent(inout) :: ipar
                double precision, intent(inout) :: x
                double precision, dimension(n), intent(inout) :: y
                double precision, dimension(:), intent(inout) :: rpar
                double precision, dimension(n), intent(inout) :: f
            end subroutine fcn
            subroutine solout(nr, xold, x, y, n, con, icomp,&
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
            end subroutine solout
        end interface
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
    end subroutine dopri5

    double precision function contd5(ii, x, con, icomp, nd)
        integer, intent(inout) :: ii
        double precision, intent(inout) :: x
        double precision, dimension(5*nd), intent(inout) :: con
        integer, dimension(nd), intent(inout) :: icomp
        integer, intent(inout) :: nd
    end function contd5
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

    real(c_double), dimension(n) :: f

    call c_f_procpointer(cfcn, fcn)
    call c_f_procpointer(csolout, solout)

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

subroutine c_dopri5(n, cfcn, x, y, xend, rtol, atol,&
        itol, csolout, iout, work, lwork, iwork,&
        liwork, rpar, ipar, idid) bind(c)
    type(c_funptr), intent(in), value :: cfcn
    type(c_funptr), intent(in), value :: csolout
    integer(c_int), intent(inout) :: n
    integer(c_int), intent(inout) :: itol
    integer(c_int), intent(inout) :: iout
    integer(c_int), intent(inout) :: lwork
    integer(c_int), intent(inout) :: liwork
    integer(c_int), dimension(:),intent(inout) :: ipar
    integer(c_int), intent(inout) :: idid
    real(c_double), intent(inout) :: xend
    real(c_double), dimension(n), intent(inout) :: rtol
    real(c_double), dimension(n), intent(inout) :: atol
    real(c_double), dimension(:), intent(inout) :: rpar
    real(c_double), dimension(lwork), intent(inout) :: work
    integer(c_int), dimension(liwork), intent(inout) :: iwork
    real(c_double), intent(inout) :: x
    real(c_double), dimension(n), intent(inout) :: y

    procedure(c_fcn), pointer :: fcn
    procedure(c_solout), pointer :: solout

    call c_f_procpointer(cfcn, fcn)
    call c_f_procpointer(csolout, solout)
    call dopri5(n, fcn, x, y, xend, rtol, atol,&
            itol, solout, iout, work, lwork, iwork,&
            liwork, rpar, ipar, idid)
end subroutine c_dopri5

function c_contd5(ii, x, con, icomp, nd) result(ret) bind(c)
    real(c_double) :: ret
    integer(c_int), intent(inout) :: ii
    real(c_double), intent(inout) :: x
    real(c_double), dimension(5*nd), intent(inout) :: con
    integer(c_int), dimension(nd), intent(inout) :: icomp
    integer(c_int), intent(inout) :: nd
    ret = contd5(ii, x, con, icomp, nd)
end function c_contd5

end module
