!!!*************************************************************
! 文件/File: eigensparse1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: eigensparse1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

SUBROUTINE egnsps1(Hmat,rowind,colptr,ntot,nnz,shift,eval,evec,nev,nconv_out,verbosity)
  implicit none

  INTEGER, intent(IN) :: ntot,nnz,nev,verbosity
  INTEGER, intent(OUT) :: nconv_out
  INTEGER :: colptr(0:ntot), rowind(nnz)
  REAL(kind=8) :: shift
  REAL(kind=8) :: Hmat(nnz)
  REAL(kind=8), TARGET :: eval(nev), evec(ntot,nev)

!     This example program is intended to illustrate the
!     simplest case of using ARPACK in considerable detail.
!     This code may be used to understand basic usage of ARPACK
!     and as a template for creating an interface to ARPACK.
!
!     This code shows how to use ARPACK to find a few eigenvalues
!     (lambda) and corresponding eigenvectors (x) for the standard
!     eigenvalue problem:
!
!                        A*x = lambda*x
!
!     where A is a n by n real nonsymmetric matrix.
!
!     The main points illustrated here are
!
!        1) How to declare sufficient memory to find NEV
!           eigenvalues of largest magnitude.  Other options
!           are available.
!
!        2) Illustration of the reverse communication interface
!           needed to utilize the top level ARPACK routine DNAUPD
!           that computes the quantities needed to construct
!           the desired eigenvalues and eigenvectors(if requested).
!
!        3) How to extract the desired eigenvalues and eigenvectors
!           using the ARPACK routine DNEUPD.
!
!     The only thing that must be supplied in order to use this
!     routine on your problem is to change the array dimensions
!     appropriately, to specify WHICH eigenvalues you want to compute
!     and to supply a matrix-vector product
!
!                         w <-  Av
!
!     in place of the call to AV( )  below.
!
!     Once usage of this routine is understood, you may wish to explore
!     the other available options to improve convergence, to solve generalized
!     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory.
!     This codes implements
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator 
!                 (Laplacian u) + rho*(du / dx)
!         on the unit square, with zero Dirichlet boundary condition.
!
!     ... OP = A  and  B = I.
!     ... Assume "call av (nx,x,y)" computes y = A*x
!     ... Use mode 1 of DNAUPD.
!
!\BeginLib
!
!\Routines called:
!     dnaupd  ARPACK reverse communication interface routine.
!     dneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     tv      Matrix vector multiplication routine that computes T*x, 
!             where T is a tridiagonal matrix.  It is used in routine
!             av.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: %Z%
! FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R%
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The maximum dimensions for all arrays are            |
!     | set here to accommodate a problem size of            |
!     | N .le. MAXN                                          |
!     |                                                      |
!     | NEV is the number of eigenvalues requested.          |
!     |     See specifications for ARPACK usage below.       |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     | You must set:                                        |
!     |                                                      |
!     | MAXN:   Maximum dimension of the A allowed.          |
!     | MAXNEV: Maximum NEV allowed.                         |
!     | MAXNCV: Maximum NCV allowed.                         |
!     %------------------------------------------------------%
!
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
  INTEGER :: iparam(11), ipntr(14)
  INTEGER, ALLOCATABLE :: ipiv(:)
  LOGICAL, ALLOCATABLE :: select(:)
  REAL(kind=8), ALLOCATABLE :: ax(:), dwrk(:,:), vwrk(:,:), workev(:), &
                            &  workl(:), workd(:), resid(:)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
  CHARACTER ::  bmat*1, which*2
  INTEGER ::  ido, n, ncv, lworkl, info, ierr, j, ishfts, maxitr, mode1, nconv
  REAL(kind=8) :: tol, sigmar, sigmai, norm
  LOGICAL :: first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
  REAL(kind=8) :: zero
  parameter         (zero = 0.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
  REAL(kind=8) :: dlapy2, dnrm2
  external          dlapy2, dnrm2, daxpy 
!
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.                    |
!     %-------------------------------------------------%
!
!      include 'debug.h'
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!

  nconv_out = 0
  n = ntot
  sigmar = shift
  sigmai = 0.d0


!
!     %-----------------------------------------------%
!     |                                               |
!     | Specifications for ARPACK usage are set       |
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |
!     |       computed.                               |
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization.                          |
!     |                                               |
!     |    3) This is a standard problem.             |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude.                      |
!     |         (indicated by which = 'LM')           |
!     |       See documentation in DNAUPD for the     |
!     |       other options SM, LR, SR, LI, SI.       |
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
  ncv = MIN(nev+20, ntot)
  bmat  = 'I'
  which = 'SR'
!
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DNAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from DNAUPD. (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DNAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  |
!     |                                                     |
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
  lworkl  = 3*ncv**2+6*ncv 
  tol    = zero 
  ido    = 0
  info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
  ALLOCATE( ax(ntot), select(ncv), dwrk(ncv,3), vwrk(ntot,ncv), &
          & workev(3*ncv), workl(lworkl), workd(3*ntot), resid(ntot) )
  ishfts = 1
  maxitr = 2000
  mode1 = 1
!
  iparam(1) = ishfts
!
  iparam(3) = maxitr
!
  iparam(7) = mode1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
   call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, vwrk, ntot, iparam, &
              &  ipntr, workd, workl, lworkl, info )

   if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
      !!call av (nx, workd(ipntr(1)), workd(ipntr(2)))
      CALL sparsemultiply(workd(ipntr(2)),workd(ipntr(1)),Hmat, &
     &                          rowind,colptr,ntot,nnz,1)
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
      go to 10
!
   endif
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
   print *, ' '
   print *, ' Error with _naupd, info = ',info
   print *, ' Check the documentation of _naupd'
   print *, ' '
!
else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine DNEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
   rvec = .true.
!
   call dneupd ( rvec, 'A', select, dwrk, dwrk(1,2), vwrk, ntot, &
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, &
     &        resid, ncv, vwrk, ntot, iparam, ipntr, workd, workl, &
     &        lworkl, ierr )
!
!        %------------------------------------------------%
!        | The real parts of the eigenvalues are returned |
!        | in the first column of the two dimensional     |
!        | array D, and the IMAGINARY part are returned   |
!        | in the second column of D.  The corresponding  |
!        | eigenvectors are returned in the first         |
!        | NCONV (= IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise,  |
!        | an orthogonal basis for the invariant subspace |
!        | corresponding to the eigenvalues in D is       |
!        | returned in V.                                 |
!        %------------------------------------------------%

  DEALLOCATE( workd )

   if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
      print *, ' '
      print *, ' Error with _neupd, info = ', ierr
      print *, ' Check the documentation of _neupd. '
      print *, ' '
!
   else
!
      first = .true.
      nconv =  iparam(5)
      do 20 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (IPARAM(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
         if (dwrk(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
            !!call av(nx, v(1,j), ax)
   CALL sparsemultiply(ax,vwrk(1,j),Hmat,rowind,colptr,ntot,nnz,1)
            call daxpy(n, -dwrk(j,1), vwrk(1,j), 1, ax, 1)
            dwrk(j,3) = dnrm2(n, ax, 1)
            dwrk(j,3) = dwrk(j,3) / abs(dwrk(j,1))
!
         else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
            !!call av(nx, v(1,j), ax)
   CALL sparsemultiply(ax,vwrk(1,j),Hmat,rowind,colptr,ntot,nnz,1)
            call daxpy(n, -dwrk(j,1), vwrk(1,j), 1, ax, 1)
            call daxpy(n, dwrk(j,2), vwrk(1,j+1), 1, ax, 1)
            dwrk(j,3) = dnrm2(n, ax, 1)
            !!call av(nx, v(1,j+1), ax)
 CALL sparsemultiply(ax,vwrk(1,j+1),Hmat,rowind,colptr,ntot,nnz,1)
            call daxpy(n, -dwrk(j,2), vwrk(1,j), 1, ax, 1)
            call daxpy(n, -dwrk(j,1), vwrk(1,j+1), 1, ax, 1)
            dwrk(j,3) = dlapy2( dwrk(j,3), dnrm2(n, ax, 1) )
            dwrk(j,3) = dwrk(j,3) / dlapy2(dwrk(j,1),dwrk(j,2))
            dwrk(j+1,3) = dwrk(j,3)
            first = .false.
         else
            first = .true.
         end if
!
 20         continue

  ALLOCATE( ipiv(nconv) )
  ipiv(:) = 0
  CALL sorteigs(ipiv,dwrk,nconv)
  DO j = 1, nconv
    eval(j) = dwrk(ipiv(j),1)
  END DO
  nconv_out = nconv

  ALLOCATE( workd(ntot) )
  workd(:) = 0.d0

  DO j = 1,nev
    evec(:,j) = vwrk(:,ipiv(j))
    IF (verbosity >= 2) THEN
	CALL sparsemultiply(workd,evec(1,j),Hmat, &
			& rowind,colptr,ntot,nnz,1)
	norm = DOT_PRODUCT(evec(:,j),workd)
	WRITE(*,*) j, ': check --> ', norm, eval(j)
    ENDIF
  END DO

  DEALLOCATE( vwrk, workd )


!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
      call dmout(6, nconv, 3, dwrk, ncv, -6, &
     &           'Ritz values (Real, Imag) and residual residuals')
   end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
   if ( info .eq. 1) then
       print *, ' '
       print *, ' Maximum number of iterations reached.'
       print *, ' '
   else if ( info .eq. 3) then
       print *, ' ' 
       print *, ' No shifts could be applied during implicit', &
     &                ' Arnoldi update, try increasing NCV.'
       print *, ' '
   end if      
!
   print *, ' '
   print *, ' _NSIMP '
   print *, ' ====== '
   print *, ' '
   print *, ' Size of the matrix is ', n
   print *, ' The number of Ritz values requested is ', nev
   print *, ' The number of Arnoldi vectors generated', &
     &            ' (NCV) is ', ncv
   print *, ' What portion of the spectrum: ', which
   print *, ' The number of converged Ritz values is ', &
     &              nconv 
   print *, ' The number of Implicit Arnoldi update', &
     &            ' iterations taken is ', iparam(3)
   print *, ' The number of OP*x is ', iparam(9)
   print *, ' The convergence criterion is ', tol
   print *, ' '
!
end if
!
!     %---------------------------%
!     | Done with program dnsimp. |
!     %---------------------------%
!
 9000 continue
!
END SUBROUTINE egnsps1
