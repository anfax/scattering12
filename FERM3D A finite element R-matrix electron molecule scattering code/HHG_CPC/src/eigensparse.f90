!!!*************************************************************
! 文件/File: eigensparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: eigensparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************



SUBROUTINE eigenfull(Hmat,Omat,rowind,colptr,nbasis,nnz,shift, &
                     & eval,evec,nev,nconv,verbosity)
  implicit none

  INTEGER, intent(IN) :: nbasis,nnz,nev,verbosity
  INTEGER, intent(OUT) :: nconv
  INTEGER :: colptr(0:nbasis), rowind(nnz)
  REAL(kind=8) :: shift
  REAL(kind=8) :: Hmat(nnz), Omat(nnz)
  REAL(kind=8), TARGET :: eval(nev), evec(nbasis,nev)

  REAL(kind=8) ::  dnrm2, dlapy2
    external          dnrm2, dlapy2

  INTEGER :: iparam(11), ipntr(11)
  INTEGER, ALLOCATABLE :: ipiv(:), rowind_save(:)
  INTEGER :: ldv, ncv, j, info, ido, lworkl, maxitr, itercount, mxitercount, ierr,nrhs
  REAL(kind=8), ALLOCATABLE :: Shftmat(:), v(:,:), d(:,:)
  REAL(kind=8), POINTER :: mx(:), resid(:), workd(:)
  REAL(kind=8), ALLOCATABLE :: workl(:)
  REAL(kind=8) :: tol, norm
  CHARACTER :: bmat*1, which*2, save_fact
  LOGICAL, ALLOCATABLE :: select(:)
  LOGICAL :: rvec 


  IF ( nev > nbasis ) THEN
    print *, ' ERROR with DSDRV4: NEV is greater than n'
    STOP 121
  END IF

  ALLOCATE( Shftmat(nnz), rowind_save(nnz) )
  rowind_save(:) = rowind(:)
  Shftmat(:) = Hmat(:) - shift*Omat(:)
  ! tol used as dummy variable, not used
  save_fact = 'Y'
  nrhs = 0
  ! LU factorize the matrix Shftmat
  CALL dSLUsolve(nbasis,nnz,nrhs,Shftmat,rowind,colptr,tol,1, &
               & info,save_fact,2)
  rowind(:) = rowind_save(:)
  Shftmat(:) = Hmat(:) - shift*Omat(:)
  DEALLOCATE( rowind_save )
  nrhs = 1

!     %-----------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The PARAMETER TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in DSAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%

  tol    = 0.d0
  ido    = 0
  info   = 0
  bmat  = 'G'
  which = 'LM'
  ldv   = nbasis
  ncv   = MIN(nev+5, nbasis)
  ALLOCATE( v(ldv,ncv) )
  v(:,:) = 0.d0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 3 of DSAUPD is used     |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%

  maxitr = 500
  mxitercount = 20*maxitr

  iparam(1) = 1
  iparam(3) = maxitr
  iparam(7) = 3


  IF (nev >= 3) THEN
    !XXX!resid => evec(:,3)
    ALLOCATE( resid(nbasis) )
  ELSE
    ALLOCATE( resid(nbasis) )
  ENDIF
  resid(1:nbasis) = 0.d0
  IF (nev >= 6) THEN
    !XXX!workd => evec(:,4)
    ALLOCATE( workd(3*nbasis) )
  ELSE
    ALLOCATE( workd(3*nbasis) )
  ENDIF
  workd(1:3*nbasis) = 0.d0


  lworkl = ncv*(ncv+8)
  ALLOCATE( workl(lworkl) )
  workl(:) = 0.d0


!     %------------------------------------------%
!     | M A I N   L O O P(Reverse communication) | 
!     %------------------------------------------%

  DO itercount = 1, mxitercount

!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by PARAMETER IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
    CALL dsaupd ( ido, bmat, nbasis, which, nev, tol, resid, &
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

      CALL sparsemultiply(workd(ipntr(2)),workd(ipntr(1)),Omat, &
                        & rowind,colptr,nbasis,nnz,1)
      CALL dSLUsolve(nbasis,nnz,nrhs,Shftmat,rowind,colptr, &
               & workd(ipntr(2)),nbasis,ierr,save_fact,0)

      IF ( ierr /= 0 ) THEN
        print*, ' '
        print*, ' ERROR with SuperLU solve in eigensparse'
        print*, ' '
        STOP 123
      END IF

!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
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
      CALL dSLUsolve(nbasis,nnz,nrhs,Shftmat,rowind,colptr, &
               & workd(ipntr(2)),nbasis,ierr,save_fact,0)

      IF ( ierr/= 0 ) THEN
        print*, ' '
        print*, ' ERROR with SuperLU solve in eigensparse'
        print*, ' '
        STOP 124
      END IF

!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%

    ELSE IF ( ido == 2) THEN

!           %---------------------------------------------%
!           |          Perform  y <--- M*x                |
!           | Need matrix vector multiplication routine   |
!           | here that takes workd(ipntr(1)) as input    |
!           | and returns the result to workd(ipntr(2)).  |
!           %---------------------------------------------%

      CALL sparsemultiply(workd(ipntr(2)),workd(ipntr(1)),Omat, &
                        & rowind,colptr,nbasis,nnz,1)

!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%

    ELSE
      EXIT
    END IF 

  END DO
  

  save_fact = 'N'
  nrhs = 0
  !! delete LU factors
  CALL dSLUsolve(nbasis,nnz,nrhs,Shftmat,rowind,colptr,tol,1, &
               & ierr,save_fact,0)
  DEALLOCATE( Shftmat )


!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%

  IF ( info .lt. 0 ) THEN

!        %--------------------------%
!        | Error message, check the |
!        | documentation in DSAUPD. |
!        %--------------------------%

    print *, ' '
    print *, ' Error with eigensparse, info = ', info
    print *, ' Check the documentation in dsaupd.'
    print *, ' '
    RETURN
  ENDIF

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%

  rvec = .true.
  ALLOCATE( d(ncv,2), select(ncv) )
  d(:,:) = 0.d0
  select(:) = .FALSE.
  CALL dseupd ( rvec, 'All', select, d, v, ldv,  &
        & shift, bmat, nbasis, which, nev, tol, &
        & resid, ncv, v, ldv, iparam, ipntr, workd, &
        & workl, lworkl, ierr )            
  DEALLOCATE( select )

  IF (nev >= 3) THEN
    !!NULLIFY( resid )
    DEALLOCATE( resid )
  ELSE
    DEALLOCATE( resid )
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
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
    print *, ' ' 
    print *, ' Error with _neupd, info = ', ierr
    print *, ' Check the documentation of _neupd. '
    print *, ' ' 

  ELSE
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


      CALL sparsemultiply(workd,v(1,j),Hmat, &
                        & rowind,colptr,nbasis,nnz,1)
      CALL sparsemultiply(workd(nbasis+1),v(1,j),Omat, &
                        & rowind,colptr,nbasis,nnz,1)
      CALL daxpy(nbasis, -d(j,1), workd(nbasis+1), 1, workd, 1)
      d(j,2) = dnrm2(nbasis, workd, 1)
      d(j,2) = d(j,2) / abs(d(j,1))
    END DO


!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%

    IF (verbosity >= 1) THEN
      CALL dmout(6, nconv, 2, d, ncv, -6, &
          & 'Ritz values and relative residuals')
    ENDIF

  END IF

  DEALLOCATE( workd )
  ALLOCATE( mx(nbasis) )
  mx(1:nbasis) = 0.d0
  DO j = 1, nev
    CALL sparsemultiply(mx,v(1,j),Omat, &
                      & rowind,colptr,nbasis,nnz,1)
    norm = DOT_PRODUCT(v(:,j),mx)
    if (abs(norm) > 1.d-5) norm = 1.d0/sqrt(norm)
    if (v(1,j) < 0.d0) norm = -norm
    v(:,j) = norm*v(:,j)
  END DO
  DEALLOCATE( mx )

  ALLOCATE( ipiv(nev) )
  ipiv(:) = 0
  CALL sorteigs(ipiv,d,nev)

  DO j = 1, nev
    eval(j) = d(ipiv(j),1)
  END DO

  ALLOCATE( workd(nbasis) )
  workd(:) = 0.d0

  DO j = 1,nev
    evec(:,j) = v(:,ipiv(j))
    IF (verbosity >= 2) THEN
      CALL sparsemultiply(workd,evec(1,j),Hmat, &
                        & rowind,colptr,nbasis,nnz,1)
      norm = DOT_PRODUCT(evec(:,j),workd)
      WRITE(*,*) j, ': check --> ', norm, eval(j)
    ENDIF
  END DO

  DEALLOCATE( d, v, workd )

  IF (verbosity >= 3) THEN
    CALL ConvergenceInfo(info, nbasis, nev, ncv, which, nconv, iparam, tol)
  ENDIF

  RETURN
END SUBROUTINE eigenfull

