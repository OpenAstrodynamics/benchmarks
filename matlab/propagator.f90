#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use dopri

    implicit none

    integer, intent(in) :: nlhs
    integer, intent(in) :: nrhs
    integer(8), dimension(:), intent(in) :: plhs
    integer(8), dimension(:), intent(in) :: prhs

    integer(8) :: mxGetPr
    integer(8) ::  mxCreateDoubleMatrix
    integer ::  mxIsNumeric
    integer(8) ::  mxGetM, mxGetN, mhandle
    integer(8) :: mxgetstring
    integer(8) :: strlen
    integer(8) :: ret
    integer(8) :: mexprintf
    integer :: mxischar


    integer :: lwork
    integer :: liwork
    integer :: iout
    integer :: idid
    integer :: itol
    integer, dimension(:), allocatable :: iwork
    integer, dimension(:), allocatable :: ipar
    integer :: lipar
    integer :: n
    double precision, dimension(:), allocatable :: rtol_
    double precision, dimension(:), allocatable :: atol_
    double precision, dimension(:), allocatable :: work
    double precision, dimension(:), allocatable :: y0

    integer, parameter :: maxhandle = 32
    character(len=maxhandle) :: handle

    if (mxischar(prhs(1)) == 1) then
        ret = mexprintf("That Works"//new_line('a'))
    end if

    ! strlen = mxGetM(prhs(1))*mxGetN(prhs(1))
    ! if (strlen > maxhandle) then
    !      call mexErrMsgIdAndTxt ('MATLAB:propagator:maxhandle', 'Max string length 32.')
    ! endif

    ! n = size(y)
    ! lwork = 11*n + 8*n + 21
    ! liwork = n + 21
    ! allocate(work(lwork))
    ! allocate(iwork(liwork))
    ! allocate(y0(n))
    ! allocate(rtol_(n))
    ! allocate(atol_(n))

    ! call dop853(n, fcn, t, y0, tend, rtol_, atol_,&
    !     itol, soldummy, iout, work, lwork, iwork,&
    !     liwork, rpar, ipar, idid)
end subroutine mexFunction
