module propagator_module

#include "fintrf.h"

implicit none

integer, parameter :: maxhandle = 32
character(len=maxhandle) :: handle

contains

! subroutine gravity(n, t, y, f, rpar, ipar)
!     integer, intent(inout) :: n
!     double precision, intent(inout) :: t
!     double precision, dimension(n), intent(inout) :: y
!     double precision, dimension(n), intent(inout) :: f
!     double precision, dimension(:), intent(inout) :: rpar
!     integer, dimension(:),intent(inout) :: ipar
!
!     double precision :: r
!
!     r = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
!     f(1) = y(4)
!     f(2) = y(5)
!     f(3) = y(6)
!
!     f(4) = -rpar(1) * y(1) / r**3
!     f(5) = -rpar(1) * y(2) / r**3
!     f(6) = -rpar(1) * y(3) / r**3
! end subroutine gravity

subroutine gravity(n, t, y, f, rpar, ipar)
    integer, intent(inout) :: n
    double precision, intent(inout) :: t
    double precision, dimension(n), intent(inout) :: y
    double precision, dimension(n), intent(inout) :: f
    double precision, dimension(:), intent(inout) :: rpar
    integer, dimension(:),intent(inout) :: ipar

    integer :: ret
    integer, parameter :: nlhs = 1
    integer, parameter :: nrhs = 2
    mwpointer, dimension(nlhs) :: plhs
    mwpointer, dimension(nrhs) :: prhs
    integer :: mexcallmatlab
    mwpointer :: mxcreatedoublematrix
    mwpointer :: mxgetpr
    mwpointer :: pr
    mwsize :: m_
    mwsize :: n_
    mwsize :: mxgetm
    mwsize :: mxgetn

    m_ = 1
    n_ = n
    prhs(1) = mxcreatedoublematrix(m_, n_, 0)
    pr = mxgetpr(prhs(1))
    call mxcopyreal8toptr(y, pr, m_*n_)
    n_ = 1
    prhs(2) = mxcreatedoublematrix(m_, n_, 0)
    pr = mxgetpr(prhs(2))
    call mxcopyreal8toptr(rpar(1), pr, m_*n_)
    ret = mexcallmatlab(nlhs, plhs, nrhs, prhs, handle)
    pr = mxgetpr(plhs(1))
    m_ = mxgetm(plhs(1))
    n_ = mxgetn(plhs(1))
    call mxcopyptrtoreal8(pr, f, m_*n_)
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

end module propagator_module
