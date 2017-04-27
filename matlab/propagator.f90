#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use dopri
    use propagator_module

    implicit none

    integer, intent(in) :: nlhs
    integer, intent(in) :: nrhs
    mwpointer, dimension(*), intent(inout) :: plhs
    mwpointer, dimension(*), intent(inout) :: prhs

    mwpointer :: mxGetPr
    mwpointer ::  mxcreatedoublematrix
    integer ::  mxIsNumeric
    mwpointer ::  mxGetM, mxGetN
    mwpointer :: mxgetstring
    mwpointer :: strlen
    mwpointer :: ret
    integer :: mxischar
    mwsize :: m
    mwsize :: n
    mwsize :: m_rpar
    mwsize :: n_rpar
    mwpointer :: out_pr
    mwpointer :: y_pr
    mwpointer :: t_pr
    mwpointer :: tend_pr
    mwpointer :: rpar_pr
    mwpointer :: locked_pr
    integer :: mexislocked

    mwsize, parameter :: one = 1

    integer :: n_
    integer :: lwork
    integer :: liwork
    integer :: iout
    integer :: idid
    integer :: itol
    integer, dimension(:), allocatable :: iwork
    integer, dimension(:), allocatable :: ipar
    integer :: lipar
    double precision, dimension(:), allocatable :: rtol
    double precision, dimension(:), allocatable :: rpar
    double precision, dimension(:), allocatable :: atol
    double precision, dimension(:), allocatable :: work
    double precision, dimension(:), allocatable :: y0
    double precision :: tend
    double precision :: t
    integer :: locked

    if (nrhs == 1) then
        locked_pr = mxgetpr(prhs(1))
        call mxcopyptrtointeger4(locked_pr, locked, one)
        if (locked == 1 .and. mexislocked() == 0) then
            call mexlock()
        elseif (locked == 0 .and. mexislocked() == 1) then
            call mexunlock()
        end if
        return
    end if

    if (nrhs /= 5) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'Five inputs required.')
    end if

    if (mxischar(prhs(1)) /= 1) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'First argument must be a function handle.')
    end if
    strlen = mxGetM(prhs(1))*mxGetN(prhs(1))
    if (strlen > maxhandle) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:maxhandle', 'Max string length 32.')
    endif
    ret = mxgetstring(prhs(1), handle, maxhandle)

    if (mxisnumeric(prhs(2)) /= 1 .and. (mxgetm(prhs(1)) == 1 .and. mxgetn(prhs(1)) == 1)) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'Second argument must be an array.')
    end if
    m = mxgetm(prhs(2))
    n = mxgetn(prhs(2))
    y_pr = mxgetpr(prhs(2))
    allocate(y0(m*n))
    call mxcopyptrtoreal8(y_pr, y0, m*n)

    if (mxisnumeric(prhs(3)) /= 1 .and. (mxgetm(prhs(1)) /= 1 .and. mxgetn(prhs(1)) /= 1)) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'Third argument must be a scalar.')
    end if
    t_pr = mxgetpr(prhs(3))
    call mxcopyptrtoreal8(t_pr, t, one)

    if (mxisnumeric(prhs(4)) /= 1 .and. (mxgetm(prhs(1)) /= 1 .and. mxgetn(prhs(1)) /= 1)) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'Fourth argument must be a scalar.')
    end if
    tend_pr = mxgetpr(prhs(4))
    call mxcopyptrtoreal8(tend_pr, t, one)

    if (mxisnumeric(prhs(5)) /= 1 .and. (mxgetm(prhs(1)) == 1 .and. mxgetn(prhs(1)) == 1)) then
         call mexErrMsgIdAndTxt ('ICATT:propagator:WrongInput', 'Fifth argument must be an array.')
    end if
    m_rpar = mxgetm(prhs(5))
    n_rpar = mxgetn(prhs(5))
    rpar_pr = mxgetpr(prhs(5))
    allocate(rpar(m_rpar*n_rpar))
    call mxcopyptrtoreal8(rpar_pr, rpar, m_rpar*n_rpar)

    n_ = n
    lwork = 11*n_ + 8*n_ + 21
    liwork = n_ + 21
    allocate(work(lwork))
    work = 0d0
    allocate(iwork(liwork))
    iwork = 0
    allocate(rtol(n_))
    rtol = 1e-6
    allocate(atol(n_))
    atol = 1e-8
    allocate(ipar(1))

    call dop853(n_, gravity, t, y0, tend, rtol, atol,&
        itol, soldummy, iout, work, lwork, iwork,&
        liwork, rpar, ipar, idid)

    if (nlhs == 1) then
        plhs(1) = mxcreatedoublematrix(m, n, 0)
        out_pr = mxgetpr(plhs(1))
        call mxcopyreal8toptr(y0, out_pr, m*n)
    else
        call mexErrMsgIdAndTxt ('ICATT:propagator:WrongOutput', 'Wrong number of output arguments.')
    end if
end subroutine mexFunction
