!!!*************************************************************
! 文件/File: dndrv4.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: dndrv4.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM dndrv4
  implicit none
!
!     Simple program to illustrate the idea of reverse communication
!     in shift-invert mode for a generalized nonsymmetric eigenvalue 
!     problem.
!
!     We implement example four of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-4
!     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
!         where A and B are derived from the finite element discretization
!         of the 1-dimensional convection-diffusion operator
!                           (d^2u / dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero Dirichlet boundary condition
!         using linear elements.
!
!     ... The shift sigma is a real number.
!
!     ... OP = inv[A-SIGMA*M]*M  and  B = M.
!
!     ... Use mode 3 of DNAUPD.
!
!\Routines called:
!     dnaupd  ARPACK reverse communication interface routine.
!     dneupd  ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors. 
!     dgttrf  LAPACK tridiagonal factorization routine.
!     dgttrs  LAPACK tridiagonal linear system solve routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     ddot    Level 1 BLAS that computes the dot product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     av      Matrix vector multiplication routine that computes A*x.
!     mv      Matrix vector multiplication routine that computes M*x.
!
!-----------------------------------------------------------------------
!
!
!     %--------------%
!     | Local Arrays |
!     %--------------%

  INTEGER :: iparam(11), ipntr(14)
  INTEGER, ALLOCATABLE :: ipiv(:)
  LOGICAL, ALLOCATABLE :: select(:)
  REAL(kind=8), ALLOCATABLE :: d(:,:), workl(:), workev(:), v(:,:), &
    & ax(:), mx(:), resid(:), workd(:), dd(:), dl(:), du(:), du2(:)

!     %---------------%
!     | Local Scalars |
!     %---------------%

  CHARACTER :: bmat*1, which*2
  INTEGER :: ido, n, nev, ncv, lworkl, info, ierr, j, ldv, &
     &                  nconv, maxitr, ishfts, mode, itercount
  REAL(kind=8) :: tol, h, s, sigmar, sigmai, s1, s2, s3
  LOGICAL :: first, rvec 
 
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%

  REAL(kind=8) :: ddot, dnrm2, dlapy2
    external          ddot, dnrm2, dlapy2, dgttrf, dgttrs

!     %------------%
!     | Parameters |
!     %------------%

  REAL(kind=8) :: one, zero, two, six, rho
    PARAMETER (one = 1.0D+0, zero = 0.0D+0, two = 2.0D+0, six = 6.0D+0)

!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues (closest   |
!     | to SIGMAR) to be approximated.  Since the          |
!     | shift-invert mode is used,  WHICH is set to 'LM'.  |
!     | The user can modify NEV, NCV, SIGMAR to solve      |
!     | problems of different sizes, and to get different  |
!     | parts of the spectrum.  However, The following     |
!     | conditions must be satisfied:                      |
!     |                     N <= MAXN,                     | 
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 2 <= NCV <= MAXNCV             | 
!     %----------------------------------------------------%

  n     = 100 
  ldv   = n
  nev   = 4 
  ncv   = 10 
  IF ( nev > n ) THEN
    print *, ' ERROR with _NDRV4: NEV is greater than n'
    STOP 121
  END IF
  ALLOCATE( ipiv(n), select(ncv) )
  ALLOCATE( d(ncv,3), workl(3*ncv*ncv+6*ncv), workev(3*ncv), v(ldv,ncv), &
    & ax(n), mx(n), resid(n), workd(3*n), dd(n), dl(n), du(n), du2(n) )
  bmat  = 'G'
  which = 'LM'
  sigmar = one 
  sigmai = zero

!     %--------------------------------------------------%
!     | Construct C = A - SIGMA*M in real arithmetic,    |
!     | and factor C in real arithmetic (using LAPACK    |
!     | SUBROUTINE dgttrf). The matrix A is chosen to be |
!     | the tridiagonal matrix derived from the standard |
!     | central difference discretization of the 1-d     |
!     | convection-diffusion operator u" + rho*u' on the |
!     | interval [0, 1] with zero Dirichlet boundary     |
!     | condition.  The matrix M is the mass matrix      |
!     | formed by using piecewise linear elements on     |
!     | [0,1].                                           |
!     %--------------------------------------------------%

  rho = 1.0D+1
  h = one / dble(n+1)
  s = rho / two

  s1 = -one/h - s - sigmar*h/six
      s2 = two/h - 4.0D+0*sigmar*h/six
      s3 = -one/h + s - sigmar*h/six

  DO j = 1, n-1
    dl(j) = s1 
    dd(j) = s2
    du(j) = s3
  END DO
  dd(n) = s2 
 
  CALL dgttrf(n, dl, dd, du, du2, ipiv, ierr)
  IF ( ierr /= 0 ) THEN
    print*, ' '
    print*, ' ERROR with _gttrf in _NDRV4.'
    print*, ' '
    STOP 122
  END IF

!     %-----------------------------------------------------%
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The PARAMETER TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in DNAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%

  lworkl = 3*ncv**2+6*ncv 
  tol    = zero
  ido    = 0
  info   = 0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 3 of DNAUPD is used     |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in DNAUPD.                          |
!     %---------------------------------------------------%

  ishfts = 1
  maxitr = 300
  mode   = 3

  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode

!     %------------------------------------------%
!     | M A I N   L O O P(Reverse communication) | 
!     %------------------------------------------%
!
  DO itercount = 1, maxitr
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take | 
!        | actions indicated by PARAMETER IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

    CALL dnaupd ( ido, bmat, n, which, nev, tol, resid, &
      & ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )

    IF (ido == -1) THEN

!           %-------------------------------------------%
!           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
!           | to force starting vector into the range   |
!           | of OP.   The user should supply his/her   |
!           | own matrix vector multiplication routine  |
!           | and a linear system solver.  The matrix   |
!           | vector multiplication routine should take |
!           | workd(ipntr(1)) as the input. The final   |
!           | result should be returned to              |
!           | workd(ipntr(2)).                          |
!           %-------------------------------------------%

      CALL mv (n, workd(ipntr(1)), workd(ipntr(2)))
      CALL dgttrs('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr) 
      IF ( ierr /= 0 ) THEN
        print*, ' '
        print*, ' ERROR with _gttrs in _NDRV4.'
        print*, ' '
        STOP 123
      END IF

!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%

    ELSE IF ( ido == 1) THEN

!           %-----------------------------------------%
!           | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
!           | M*x has been saved in workd(ipntr(3)).  |
!           | The user only need the linear system    |
!           | solver here that takes workd(ipntr(3))  |
!           | as input, and returns the result to     |
!           | workd(ipntr(2)).                        |
!           %-----------------------------------------%

      CALL dcopy( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
      CALL dgttrs ('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr)
      IF ( ierr /= 0 ) THEN
        print*, ' '
        print*, ' ERROR with _gttrs in _NDRV4.'
        print*, ' '
        STOP 124
      END IF

!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%

    ELSE IF ( ido == 2) THEN

!           %---------------------------------------------%
!           |          Perform  y <--- M*x                |
!           | Need matrix vector multiplication routine   |
!           | here that takes workd(ipntr(1)) as input    |
!           | and returns the result to workd(ipntr(2)).  |
!           %---------------------------------------------%

      CALL mv (n, workd(ipntr(1)), workd(ipntr(2)))

!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%

    ELSE

      EXIT

    END IF 

  END DO


!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%

    IF ( info .lt. 0 ) THEN

!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%

      print *, ' '
      print *, ' Error with _naupd, info = ', info
      print *, ' Check the documentation in _naupd.'
      print *, ' '

    ELSE 

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%

      rvec = .true.
      CALL dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,  &
        & sigmar, sigmai, workev, bmat, n, which, nev, tol, &
        & resid, ncv, v, ldv, iparam, ipntr, workd, &
        & workl, lworkl, ierr )            

!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the IMAGINARY part is returned   |
!        | in the second column of D.  The corresponding |
!        | eigenvectors are returned in the first NEV    |
!        | columns of the two dimensional array V if     |
!        | requested.  Otherwise, an orthogonal basis    |
!        | for the invariant subspace corresponding to   |
!        | the eigenvalues in D is returned in V.        |
!        %-----------------------------------------------%

      IF ( ierr /= 0) THEN
 
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DNEUPD. |
!            %------------------------------------%
!
        print *, ' ' 
        print *, ' Error with _neupd, info = ', ierr
        print *, ' Check the documentation of _neupd. '
        print *, ' ' 
      ELSE
        first = .true.
        nconv =  iparam(5)
        DO j=1, nconv

!                %---------------------------%
!                | Compute the residual norm |
!                |                           |
!                |   ||  A*x - lambda*x ||   |
!                |                           |
!                | for the NCONV accurately  |
!                | computed eigenvalues and  |
!                | eigenvectors.  (iparam(5) |
!                | indicates how many are    |
!                | accurate to the requested |
!                | tolerance)                |
!                %---------------------------%
!
          IF (d(j,2) == zero)  THEN
!
!                    %--------------------%
!                    | Ritz value is real |
!                    %--------------------%

            CALL av(n, v(1,j), ax, rho)
            CALL mv(n, v(1,j), mx)
            CALL daxpy(n, -d(j,1), mx, 1, ax, 1)
            d(j,3) = dnrm2(n, ax, 1)
            d(j,3) = d(j,3) / abs(d(j,1))

          ELSE IF (first) THEN

!                    %------------------------%
!                    | Ritz value is complex. |
!                    | Residual of one Ritz   |
!                    | value of the conjugate |
!                    | pair is computed.      | 
!                    %------------------------%
!     
            CALL av(n, v(1,j), ax, rho)
            CALL mv(n, v(1,j), mx)
            CALL daxpy(n, -d(j,1), mx, 1, ax, 1)
            CALL mv(n, v(1,j+1), mx)
            CALL daxpy(n, d(j,2), mx, 1, ax, 1)
            d(j,3) = dnrm2(n, ax, 1)
            CALL av(n, v(1,j+1), ax, rho)
            CALL mv(n, v(1,j+1), mx)
            CALL daxpy(n, -d(j,1), mx, 1, ax, 1)
            CALL mv(n, v(1,j), mx)
            CALL daxpy(n, -d(j,2), mx, 1, ax, 1)
            d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
            d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
            d(j+1,3) = d(j,3)
            first = .false.
          ELSE
            first = .true.
          END IF
        END DO

!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%

        CALL dmout(6, nconv, 3, d, ncv, -6, &
          & 'Ritz values (Real,Imag) and relative residuals')

      END IF

!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%

      IF ( info == 1) THEN
        print *, ' '
        print *, ' Maximum number of iterations reached.'
        print *, ' '
      ELSE IF ( info == 3) THEN
        print *, ' ' 
        print *, ' No shifts could be applied during implicit', &
          & ' Arnoldi update, try increasing NCV.'
        print *, ' '
      END IF      

      print *, ' '
      print *, ' _NDRV4 '
      print *, ' ====== '
      print *, ' '
      print *, ' Size of the matrix is ', n
      print *, ' The number of Ritz values requested is ', nev
      print *, ' The number of Arnoldi vectors generated', &
        & ' (NCV) is ', ncv
      print *, ' What portion of the spectrum: ', which
      print *, ' The number of converged Ritz values is ', &
        & nconv 
      print *, ' The number of Implicit Arnoldi update', &
        & ' iterations taken is ', iparam(3)
      print *, ' The number of OP*x is ', iparam(9)
      print *, ' The convergence criterion is ', tol
      print *, ' '

    END IF

    STOP
  END PROGRAM dndrv4
 

!==========================================================================

!     matrix vector multiplication SUBROUTINE

SUBROUTINE mv (n, v, w)
  implicit none
  INTEGER :: n, j
  REAL(kind=8) :: v(n), w(n), one, four, six, h 
      PARAMETER (one = 1.0D+0, four = 4.0D+0, six = 6.0D+0)

!     Compute the matrix vector multiplication y<---M*x
!     where M is mass matrix formed by using piecewise linear elements 
!     on [0,1].
 
  w(1) =  ( four*v(1) + one*v(2) ) / six
  DO j = 2,n-1
    w(j) = ( one*v(j-1) + four*v(j) + one*v(j+1) ) / six
  END DO
  w(n) =  ( one*v(n-1) + four*v(n) ) / six

  h = one / dble(n+1)
  CALL dscal(n, h, w, 1)
  RETURN
END SUBROUTINE mv


!------------------------------------------------------------------
SUBROUTINE av (n, v, w, rho)
  implicit none
  INTEGER :: n, j
  REAL(kind=8) :: v(n), w(n), one, two, dd, dl, du, s, h, rho 
      PARAMETER         (one = 1.0D+0, two = 2.0D+0)

!     Compute the matrix vector multiplication y<---A*x
!     where A is obtained from the finite element discretization of the
!     1-dimensional convection diffusion operator
!                     d^u/dx^2 + rho*(du/dx)
!     on the interval [0,1] with zero Dirichlet boundary condition
!     using linear elements.
!     This routine is only used in residual calculation.

  h = one / dble(n+1)
  s = rho / two
  dd = two / h 
  dl = -one/h - s
  du = -one/h + s

  w(1) =  dd*v(1) + du*v(2)
  DO j = 2,n-1
    w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
  END DO
  w(n) =  dl*v(n-1) + dd*v(n)
  RETURN
END SUBROUTINE av

