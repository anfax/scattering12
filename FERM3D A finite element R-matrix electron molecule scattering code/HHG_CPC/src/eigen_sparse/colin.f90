!!!*************************************************************
! 文件/File: colin.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: colin.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM colin
  USE FiniteElements
  implicit none

  REAL(kind=8), EXTERNAL :: efield

  TYPE(FEBnodes) :: x
  REAL(kind=8), ALLOCATABLE :: Hmat(:), Omat(:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:), cmask(:)
  REAL(kind=8), ALLOCATABLE :: psico(:,:), psi0(:,:), phi(:), pulse(:,:)
  REAL(kind=8), ALLOCATABLE :: O1d(:,:), H1d(:,:)
  REAL(kind=8), ALLOCATABLE :: f(:,:), r(:), savereal(:,:), saveimag(:,:)
  INTEGER :: nnz, nb2, nb1, nsteps, nbs, icbs, nsave
  INTEGER :: nev, n, j, l, nconv, col, nrhs, info, ngrd, nrghtmsk
  INTEGER, ALLOCATABLE :: rowind(:), colptr(:)
  REAL(kind=8) :: xmax, sigmar, norm, lamdamax, lamdamin, pi, deltat, dt, spect(2), alpha, beta
  CHARACTER :: save_fact, pot1_type*10, pot2_type*10, flg

  nsave = 50     ! number of snapshots to save
  deltat = 1.d0  ! time between each snapshot
  nsteps = 3     ! number of propagation steps for each snapshot

  pi = 4.d0*atan(1.d0)

  CALL CreateNodes(x,'colin.nds')
  CALL PrintNodes(x,6)

  lamdamin = abs(x%nodes(2) - x%nodes(1))
  DO j = 3,x%numnds
    lamdamax = abs(x%nodes(j) - x%nodes(j-1))
    if (lamdamax < lamdamin) lamdamin = lamdamax
  END DO
  lamdamax = 14.d0*pi*pi/(lamdamin**2)
  lamdamin = -2.d0
  spect(1) = lamdamin
  spect(2) = lamdamax

  nb1 = x%numbasis
  nb2 = nb1*nb1
  IF (x%numnds <= 2) THEN
    nnz = nb2
  ELSE
    nnz = x%numleft*(x%numleft+6) + x%numright*(x%numright+6) + &
	& 9*(3*x%numnds-8)
  ENDIF
  nnz = nnz*nnz

  ALLOCATE( O1d(-5:5,nb1), H1d(-5:5,nb1) )
  pot1_type = 'coulomb'
  alpha = -1.d0
  beta = 0.d0
  CALL Vmtrx(H1d,pot1_type,alpha,beta,x)
  CALL KEmtrx(O1d,x)
  H1d(:,:) = H1d(:,:) + O1d(:,:)
  CALL Omtrx(O1d,x)

  nev = 6
  ALLOCATE( phi(nb1), pulse(nb1,2), eval(nev), evec(nb1,nev) )
  sigmar = -.5d0
  CALL eigenfew(H1d,O1d,nb1,sigmar,eval,evec,nev,nconv,0)
  nbs = 0
  DO j = 1,nconv
    WRITE(*,*) j, eval(j)
    if (eval(j) < 0.d0) nbs = nbs + 1
  END DO
  icbs = 1
  phi(:) = evec(:,icbs)
  WRITE(*,*) 'initial state energy = ', eval(icbs)

  CALL wvpckt(pulse,x)

  ngrd = 175
  ALLOCATE( r(ngrd), f(ngrd,3) )
  DO j = 1,ngrd
    r(j) = j*0.2d0
  END DO
  CALL funeval1(f(1,1),r,ngrd,pulse(1,1),x)
  CALL funeval1(f(1,2),r,ngrd,pulse(1,2),x)
  CALL funeval1(f(1,3),r,ngrd,phi,x)
  OPEN(UNIT=8, FILE='test.dat')
  DO j = 1,ngrd
    WRITE(8,100) r(j), f(j,1), f(j,2), f(j,3)
  END DO
  CLOSE(8)
100 FORMAT(4E15.7)
  DEALLOCATE( f )

  ALLOCATE( psico(nb2,2) )

  CALL makeinitcon( psico(1,1), phi, pulse(1,1), nb1, 1)
  CALL makeinitcon( psico(1,2), phi, pulse(1,2), nb1, 1)

  WRITE(*,*) 'guessed number of nonzeros is ', nnz
  ALLOCATE( Hmat(nnz), rowind(nnz), colptr(0:nb2), Omat(nnz) )
  CALL creatematrices(Hmat, Omat, rowind, colptr, nnz, H1d, O1d, x)
  WRITE(*,*) 'actual number of nonzeros is ', nnz

  ALLOCATE( savereal(nb2,0:nsave), saveimag(nb2,0:nsave), cmask(nb2) )
  DO nrghtmsk = 1,x%numnds
    if ( x%nodes(x%numnds-nrghtmsk) < (x%nodes(x%numnds)-30.d0) ) EXIT
  END DO
  CALL coefmask(savereal,x,0,nrghtmsk)
  savereal(1:nb1,0) = savereal(1:nb1,0)**2
  CALL makeinitcon( cmask, savereal, savereal, nb1, 0)

  CALL sparsemultiply(savereal(1,0),psico,Omat,rowind,colptr,nb2,nnz,2)
  norm = DOT_PRODUCT(psico(:,1),savereal(:,0)) + DOT_PRODUCT(psico(:,2),savereal(:,1))
  WRITE(*,*) 'norm of initial state is ', norm

  !! save initial state
  savereal(:,0) = psico(:,1)
  saveimag(:,0) = psico(:,2)
  WRITE(*,*) 't = ', 0*deltat
  dt = deltat/nsteps

  flg = 'Y'
  DO j = 1,nsave
    DO l = 1,nsteps
      psico(:,1) = cmask(:)*psico(:,1)
      psico(:,2) = cmask(:)*psico(:,2)
      if ( (j==nsave).AND.(l==nsteps) ) flg = 'N'
      CALL chebysparse(psico,dt,Hmat,Omat,rowind,colptr,nb2,nnz,spect,flg)
    END DO
    WRITE(*,*) 't = ', j*deltat
    savereal(:,j) = psico(:,1)
    saveimag(:,j) = psico(:,2)
  END DO

  OPEN(UNIT=8, FILE='colin.bin', FORM='UNFORMATTED')
  WRITE(8) nb1, nsave, savereal, saveimag
  CLOSE(8)
  WRITE(*,*) 'saved data to colin.bin'

  DO j = 0,nsave
    CALL projectout(savereal(1,j),evec,O1d,nb1,nbs)
    CALL projectout(saveimag(1,j),evec,O1d,nb1,nbs)
  END DO

  DEALLOCATE( evec )
  ALLOCATE( f(ngrd*ngrd,2), evec(ngrd*ngrd,0:nsave) )
  DO j = 0,nsave
    CALL funeval2(f(1,1),r,ngrd,savereal(1,j),x)
    CALL funeval2(f(1,2),r,ngrd,saveimag(1,j),x)
    evec(:,j) = f(:,1)*f(:,1) + f(:,2)*f(:,2)
  END DO

  STOP

  OPEN(UNIT=8, FILE='t2d.dat')
  DO j = 1,ngrd*ngrd
    WRITE(8,200) evec(j,:)
  END DO
  CLOSE(8)
200 FORMAT(100E11.3)


  STOP
END PROGRAM colin


SUBROUTINE projectout(psi,phi,O1d,nb1,nbs)
  implicit none

  INTEGER :: nb1, nbs
  REAL(kind=8) :: psi(nb1,nb1)
  REAL(kind=8) :: O1d(-5:5,nb1), phi(nb1,nbs)
  INTEGER :: j, k, mu
  REAL(kind=8), ALLOCATABLE :: p1(:,:), p2(:,:), p12(:), work(:)

  ALLOCATE( p1(nb1,nbs), p2(nb1,nbs), p12(nbs), work(nb1) )

  DO j = 1,nb1
    p1(:,1) = psi(j,:)
    CALL Av_FEband(nb1,O1d,p1(1,1),work,1)
    DO mu = 1,nbs
      p2(j,mu) = DOT_PRODUCT(phi(:,mu),work(:))
    END DO
  END DO

  DO mu = 1,nbs
    CALL Av_FEband(nb1,O1d,p2(1,mu),work,1)
    p12(mu) = DOT_PRODUCT(phi(:,mu),work(:))
  END DO

  DO j = 1,nb1
    CALL Av_FEband(nb1,O1d,psi(1,j),work,1)
    DO mu = 1,nbs
      p1(j,mu) = DOT_PRODUCT(phi(:,mu),work(:))
    END DO
  END DO

  DO j = 1,nb1
    DO k = 1,nb1
      DO mu = 1,nbs
        psi(j,k) = psi(j,k) -  phi(j,mu)*p1(k,mu) - p2(j,mu)*phi(k,mu) + p12(mu)*phi(j,mu)*phi(k,mu)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE projectout


SUBROUTINE creatematrices(Hsparse, Osparse, rowind, colptr, nnzinout, H1d, O1d, x)
  USE FiniteElements
  USE Potentials
  implicit none

  INTEGER :: nnzinout
  TYPE(FEBnodes) :: x
  INTEGER :: rowind(nnzinout), colptr(0:*)
  REAL(kind=8) :: Hsparse(nnzinout), Osparse(nnzinout)
  REAL(kind=8) :: H1d(-5:5,x%numbasis),  O1d(-5:5,x%numbasis)
  INTEGER :: nbasis, nnz, nb2
  REAL(kind=8), ALLOCATABLE :: H2d(:,:), O2d(:,:)
  REAL(kind=8) :: alpha, beta
  CHARACTER :: pot1_type*10, pot2_type*10

  alpha = -1.d0
  beta = 0.d0
  pot2_type = 'colinear'

  nbasis = x%numbasis
  nb2 = nbasis*nbasis

  ALLOCATE( H2d(-60:60,nb2), O2d(-60:60,nb2) )
  H2d(:,:) = 0.d0
  O2d(:,:) = 0.d0

  CALL ham2d(H2d,H1d,O1d,H1d,O1d,nbasis,nbasis)
  !! O2d temporarily holds the 2-electron potential
  CALL V2dmtrx(O2d,pot2_type,x,x)
  H2d(:,:) = H2d(:,:) + O2d(:,:)

  O2d(:,:) = 0.d0
  CALL directsum(O2d,O1d,O1d,nbasis,nbasis)

  CALL sparse2d(Hsparse,rowind,colptr,nb2,nnz,H2d,nbasis,nbasis)

  CALL sparsegrab2d(Osparse,rowind,colptr,nb2,nnz,O2d,nbasis,nbasis)

  nnzinout = nnz
  RETURN
END SUBROUTINE creatematrices


SUBROUTINE makeinitcon( psi0, f1, f2, nb1, sym)
  implicit none

  INTEGER :: nb1, sym
  REAL(kind=8) :: psi0(nb1,nb1), f1(nb1), f2(nb1)
  INTEGER :: j, k
  REAL(kind=8) :: norm

  IF ( sym == 0 ) THEN
    DO k = 1,nb1
      DO j = 1,nb1
        psi0(j,k) = f1(j)*f2(k)
      END DO
    END DO
  ELSE
    norm = 1.d0 /sqrt(2.d0)
    DO k = 1,nb1
      DO j = 1,nb1
        psi0(j,k) = norm*(f1(j)*f2(k) + sym*f2(j)*f1(k))
      END DO
    END DO
  ENDIF

  RETURN
END SUBROUTINE makeinitcon
