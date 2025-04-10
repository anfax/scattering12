!!!*************************************************************
! 文件/File: eigen2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: eigen2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************



SUBROUTINE eigen2d(Hmat,Omat,nb1,nb2,sigmar,eval,evec,nev,nconv,verbosity)
  implicit none

  INTEGER :: nb1, nb2, nev, nconv, verbosity
  REAL(kind=8) :: sigmar
  REAL(kind=8) :: Hmat(-60:60,nb1*nb2), Omat(-60:60,nb1*nb2)
  REAL(kind=8), TARGET :: eval(nev), evec(nb1*nb2,nev)

  REAL(kind=8) ::  dnrm2, dlapy2
    external          dnrm2, dlapy2

  INTEGER :: iparam(11), ipntr(14)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER :: ldv, ncv, ldab, j, info, ido, nrhs, lworkl, ishfts, maxitr, mode, itercount, ierr, nbasis, kl, ku
  REAL(kind=8), ALLOCATABLE :: Shftmat(:,:), v(:,:), d(:,:)
  REAL(kind=8), POINTER :: ax(:), mx(:), resid(:), workd(:)
  REAL(kind=8), ALLOCATABLE :: workl(:), workev(:)
  REAL(kind=8) :: sigmai, tol, norm
  CHARACTER :: bmat*1, which*2
  LOGICAL, ALLOCATABLE :: select(:)
  LOGICAL :: first, rvec 

  sigmai = 0.d0
  IF ( nev > nb1*nb2 ) THEN
    print *, ' ERROR with _NDRV4: NEV is greater than n'
    STOP 121
  END IF

  kl = nb1*min(5,nb2-1) + min(5,nb1-1)
  ku = kl
  nbasis = nb1*nb2
  ldab = 2*kl + ku + 1
  ALLOCATE( Shftmat(ldab,nbasis) )
  CALL shftmat2d(Shftmat,Hmat,Omat,sigmar,nb1,nb2,ldab)

  ALLOCATE( ipiv(nbasis) )
  CALL dgbtrf(nbasis,nbasis,kl,ku,Shftmat,ldab,ipiv,info)

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

  tol    = 0.d0
  ido    = 0
  info   = 0
  nrhs = 1
  bmat  = 'G'
  which = 'LM'
  ldv   = nbasis
  ncv   = MIN(nev+5, nbasis)
  ALLOCATE( v(ldv,ncv) )

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

  ax => evec(:,1)
  IF (nev >= 2) THEN
    mx => evec(:,2)
  ELSE
    ALLOCATE( mx(nbasis) )
  ENDIF
  IF (nev >= 3) THEN
    resid => evec(:,3)
  ELSE
    ALLOCATE( resid(nbasis) )
  ENDIF
  IF (nev >= 6) THEN
    workd => evec(:,4)
  ELSE
    ALLOCATE( workd(3*nbasis) )
  ENDIF

  lworkl = 3*ncv**2+6*ncv 
  ALLOCATE( workl(lworkl) )

  DO itercount = 1, maxitr

!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take | 
!        | actions indicated by PARAMETER IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

    CALL dnaupd ( ido, bmat, nbasis, which, nev, tol, resid, &
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

      !CALL Av_FEband(nbasis, Omat, workd(ipntr(1)), workd(ipntr(2)), 1)
      CALL Av_FE2band(Omat,nb1,nb2, workd(ipntr(1)), workd(ipntr(2)), 1)
      CALL dgbtrs('N',nbasis,kl,ku,nrhs,Shftmat,ldab,ipiv,workd(ipntr(2)),nbasis,ierr)
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

      CALL dcopy( nbasis, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
      CALL dgbtrs('N',nbasis,kl,ku,nrhs,Shftmat,ldab,ipiv,workd(ipntr(2)),nbasis,ierr)
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

      !CALL Av_FEband(nbasis, Omat, workd(ipntr(1)), workd(ipntr(2)), 1)
      CALL Av_FE2band(Omat,nb1,nb2, workd(ipntr(1)), workd(ipntr(2)), 1)

!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%

    ELSE

      EXIT

    END IF 

  END DO

  DEALLOCATE( Shftmat, ipiv )


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
    RETURN
  ENDIF

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
  ALLOCATE( d(ncv,3) )

  ALLOCATE( workev(3*ncv), select(ncv) )
  CALL dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,  &
        & sigmar, sigmai, workev, bmat, nbasis, which, nev, tol, &
        & resid, ncv, v, ldv, iparam, ipntr, workd, &
        & workl, lworkl, ierr )            
  DEALLOCATE( workev, select )

  IF (nev >= 3) THEN
    NULLIFY( resid )
  ELSE
    DEALLOCATE( resid )
  ENDIF
  IF (nev >= 6) THEN
    NULLIFY( workd )
  ELSE
    DEALLOCATE( workd )
  ENDIF
  DEALLOCATE( workl )

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

      IF (d(j,2) == 0.d0)  THEN

!                    %--------------------%
!                    | Ritz value is real |
!                    %--------------------%

        !CALL Av_FEband(nbasis, Hmat, v(1,j), ax, 1)
        CALL Av_FE2band(Hmat,nb1,nb2, v(1,j), ax, 1)
        !CALL Av_FEband(nbasis, Omat, v(1,j), mx, 1)
        CALL Av_FE2band(Omat,nb1,nb2, v(1,j), mx, 1)
        CALL daxpy(nbasis, -d(j,1), mx, 1, ax, 1)
        d(j,3) = dnrm2(nbasis, ax, 1)
        d(j,3) = d(j,3) / abs(d(j,1))

      ELSE IF (first) THEN

!                    %------------------------%
!                    | Ritz value is complex. |
!                    | Residual of one Ritz   |
!                    | value of the conjugate |
!                    | pair is computed.      | 
!                    %------------------------%

        !CALL Av_FEband(nbasis, Hmat, v(1,j), ax, 1)
        CALL Av_FE2band(Hmat,nb1,nb2, v(1,j), ax, 1)
        !CALL Av_FEband(nbasis, Omat, v(1,j), mx, 1)
        CALL Av_FE2band(Omat,nb1,nb2, v(1,j), mx, 1)
        CALL daxpy(nbasis, -d(j,1), mx, 1, ax, 1)
        !CALL Av_FEband(nbasis, Omat, v(1,j+1), mx, 1)
        CALL Av_FE2band(Omat,nb1,nb2, v(1,j+1), mx, 1)
        CALL daxpy(nbasis, d(j,2), mx, 1, ax, 1)
        d(j,3) = dnrm2(nbasis, ax, 1)
        !CALL Av_FEband(nbasis, Hmat, v(1,j+1), ax, 1)
        CALL Av_FE2band(Hmat,nb1,nb2, v(1,j+1), ax, 1)
        !CALL Av_FEband(nbasis, Omat, v(1,j+1), mx, 1)
        CALL Av_FE2band(Omat,nb1,nb2, v(1,j+1), mx, 1)
        CALL daxpy(nbasis, -d(j,1), mx, 1, ax, 1)
        !CALL Av_FEband(nbasis, Omat, v(1,j), mx, 1)
        CALL Av_FE2band(Omat,nb1,nb2, v(1,j), mx, 1)
        CALL daxpy(nbasis, -d(j,2), mx, 1, ax, 1)
        d(j,3) = dlapy2( d(j,3), dnrm2(nbasis, ax, 1) )
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

    IF (verbosity >= 1) THEN
      CALL dmout(6, nconv, 3, d, ncv, -6, &
          & 'Ritz values (Real,Imag) and relative residuals')
    ENDIF

  END IF

  DO j = 1, nev
    !CALL Av_FEband(nbasis, Omat, v(:,j), ax, 1)
    CALL Av_FE2band(Omat,nb1,nb2, v(1,j), ax, 1)
    norm = DOT_PRODUCT(v(:,j),ax)
    norm = 1.d0/sqrt(norm)
    if (v(1,j) < 0.d0) norm = -norm
    v(:,j) = norm*v(:,j)
  END DO

  ALLOCATE( ipiv(nev) )
  CALL sorteigs(ipiv,d,nev)

! Note, this returns only the real part of the eigenvalue!!!
  DO j = 1, nev
    eval(j) = d(ipiv(j),1)
  END DO

  NULLIFY( ax )
  IF (nev >= 2) THEN
    NULLIFY( mx )
  ELSE
    DEALLOCATE( mx )
  ENDIF
  ALLOCATE( workd(nbasis) )

  DO j = 1,nev
    evec(:,j) = v(:,ipiv(j))
    IF (verbosity >= 2) THEN
      !CALL Av_FEband(nbasis, Hmat, evec(:,j), workd, 1)
      CALL Av_FE2band(Hmat,nb1,nb2, evec(:,j), workd, 1)
      norm = DOT_PRODUCT(evec(:,j),workd)
      WRITE(*,*) j, ': check --> ', norm, eval(j)
    ENDIF
  END DO

  DEALLOCATE( d, v )

  IF (verbosity >= 3) THEN
    CALL ConvergenceInfo(info, nbasis, nev, ncv, which, nconv, iparam, tol)
  ENDIF

  RETURN
END SUBROUTINE eigen2d


SUBROUTINE shftmat2d(shftmat,hmat,omat,a,nb1,nb2,ldab)
  implicit none

  INTEGER :: nb1, nb2, ldab
  REAL(kind=8) :: shftmat(ldab,*), hmat(-60:60,nb1*nb2), omat(-60:60,nb1*nb2)
  REAL(kind=8) :: a
  INTEGER :: col, blk, id, jd, kl, ku, strd1, strd2, centroid, j, d1, d2

  d1 = min(5,nb1-1)
  d2 = min(5,nb2-1)
  kl = d2*nb1 + d1
  ku = kl
  IF ( ldab < (2*kl + ku + 1) ) THEN
    WRITE(*,*) 'fatal error in shftmat2d: ldab (=', ldab, ') is too small'
    STOP 181
  ENDIF
  centroid = kl + ku + 1
  strd1 = 2*d1 + 1
  strd2 = MAX(nb1,strd1)

  DO col = 1, nb1*nb2
    DO blk = -d2,d2
      id = blk*strd2 + centroid
      jd = blk*strd1
      shftmat(id-d1:id+d1,col) = hmat(jd-d1:jd+d1,col) - a*omat(jd-d1:jd+d1,col)
    END DO
  END DO

  RETURN
END SUBROUTINE shftmat2d

