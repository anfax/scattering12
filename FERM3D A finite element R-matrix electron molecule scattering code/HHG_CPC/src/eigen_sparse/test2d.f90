!!!*************************************************************
! 文件/File: test2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: test2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM test2d
  USE FiniteElements
  USE potentials
  implicit none

  TYPE(FEBnodes) :: x
  REAL(kind=8), ALLOCATABLE :: Hmat(:,:), Omat(:,:), Imat(:,:)
  REAL(kind=8), ALLOCATABLE :: h2d(:,:), o2d(:,:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:), cmask(:)
  INTEGER :: nmnds, nbasis, nsave, nprop
  INTEGER :: j, nconv, ldab, info, ierr, coarsesteps, nrhs
  INTEGER, ALLOCATABLE :: ipiv(:)
  REAL(kind=8) :: xmax, sigmar, t0, E0, savetime, t, dt, norm, imag, Amp, error

  t0 = 0.d0
  E0 = 1.d-3
  nsave = 570
  savetime = 1.d0
  nprop = 100

  CALL CreateNodes(x,'nodes.inp')
  CALL PrintNodes(x,6)
  nmnds = x%numnds

  nbasis = 3*(nmnds-2) + 4
  ALLOCATE( Omat(-5:5,nbasis), Hmat(-5:5,nbasis) )
  CALL Omtrx(Omat,x)
  CALL KEmtrx(Hmat,x)


  ALLOCATE( o2d(-60:60,nbasis*nbasis), h2d(-60:60,nbasis*nbasis) )
  o2d(:,:) = 0.d0
  h2d(:,:) = 0.d0
  ALLOCATE( Imat(-5:5,nbasis) )
  Imat(:,:) = 0.d0
  !Omat(:,:) = 0.d0
  !Omat(-5,6:nbasis) = -5.d0
  !Imat(-5,6:nbasis) = 1.d0
  !Omat(-4,5:nbasis) = -4.d0
  !Imat(-4,5:nbasis) = 1.d0
  !Omat(-3,4:nbasis) = -3.d0
  !Imat(-3,4:nbasis) = 1.d0
  !Omat(-2,3:nbasis) = -2.d0
  !Imat(-2,3:nbasis) = 1.d0
  !Omat(-1,2:nbasis) = -1.d0
  !Imat(-1,2:nbasis) = 1.d0
  !Omat(0,1:nbasis) = 1.d0
  Imat(0,1:nbasis) = 1.d0
  !Omat(1,1:nbasis-1) = 1.d0
  !Imat(1,1:nbasis-1) = 1.d0
  !Omat(2,1:nbasis-2) = 2.d0
  !Imat(2,1:nbasis-2) = 1.d0
  !Omat(3,1:nbasis-3) = 3.d0
  !Imat(3,1:nbasis-3) = 1.d0
  !Omat(4,1:nbasis-4) = 4.d0
  !Imat(4,1:nbasis-4) = 1.d0
  !Omat(5,1:nbasis-5) = 5.d0
  !Imat(5,1:nbasis-5) = 1.d0
  DO j = -5,5
  !  WRITE(*,101) j, Omat(j,:)
  END DO
101 FORMAT(I4,20E13.3)
  CALL directsum(o2d,Omat,Omat,nbasis,nbasis)
  CALL ham2d(h2d,Hmat,Omat,Hmat,Omat,nbasis,nbasis)

  h2d = 0.5d0*h2d
  !CALL dumpbandmatrix(o2d,121,nbasis*nbasis,60,60)

  ALLOCATE( evec(nbasis*nbasis,3) )
  CALL dumbshit(h2d,o2d,nbasis,nbasis,evec)
  evec(:,2) = 0.d0
  !evec(1:,1) = 1.d0
  !evec(2:nbasis*nbasis,1) = 0.d0
  !CALL matvec2d(evec(1,3),Omat,evec(1,1),nbasis)
  CALL Av_FE2band(h2d,nbasis,nbasis,evec(1,1),evec(1,2),1)
  DO j = 1,nbasis*nbasis
  !  WRITE(*,*) j, evec(j,:)
  END DO

  evec(:,1) = ABS(evec(:,2) - evec(:,3))
  error = SUM(evec(:,1))/(nbasis*nbasis)
  WRITE(*,*) 'error = ', error

  STOP
END PROGRAM test2d


SUBROUTINE dumbshit(Hmat,Omat,nb1,nb2,evec)
  implicit none

  INTEGER :: nb1, nb2
  REAL(kind=8) :: Hmat(-60:60,nb1*nb2)
  REAL(kind=8) :: Omat(-60:60,nb1*nb2)
  REAL(kind=8) :: evec(nb1*nb2,3)

  INTEGER :: iparam(11), ipntr(14)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER :: ldv, ncv, ldab, j, info, ido, nrhs, lworkl, ishfts, maxitr, mode, itercount, ierr, nbasis, kl, ku
  REAL(kind=8), ALLOCATABLE :: Shftmat(:,:)
  REAL(kind=8) :: sigmar, tol, norm
  CHARACTER :: bmat*1, which*2
  LOGICAL, ALLOCATABLE :: select(:)
  LOGICAL :: first, rvec 

  kl = nb1*min(5,nb2-1) + min(5,nb1-1)
  ku = kl
  nbasis = nb1*nb2
  ldab = 2*kl + ku + 1
  ALLOCATE( Shftmat(ldab,nbasis) )
  sigmar = 1.d0
  CALL shftmat2dumb(Shftmat,Hmat,Omat,sigmar,nb1,nb2,ldab)
  CALL dumpbandmatrix(Shftmat,ldab,nbasis,kl,ku)
  DO j = 1,nbasis
    evec(j,1) = sin(.01d0*j*j)
    evec(j,2) = evec(j,1)
    evec(j,3) = evec(j,1)
  END DO

  ALLOCATE( ipiv(nbasis) )
  CALL dgbtrf(nbasis,nbasis,kl,ku,Shftmat,ldab,ipiv,info)
  nrhs = 1
  CALL dgbtrs('N',nbasis,kl,ku,nrhs,Shftmat,ldab,ipiv,evec,nbasis,ierr)

  RETURN
END SUBROUTINE dumbshit


SUBROUTINE matvec2d(b,A1d,x,n)
  implicit none

  INTEGER :: n, nrhs
  REAL(kind=8) :: b(n,n), x(n,n), A1d(-5:5,n)
  INTEGER :: i, j
  REAL(kind=8), ALLOCATABLE :: temp(:,:)


  nrhs = n
  ALLOCATE( temp(n,n) )
  CALL Av_FEband(n, A1d, x, temp, nrhs)

  DO j = 1,n
    DO i = 1,n
      b(j,i) = temp(i,j)
    END DO
  END DO
  CALL Av_FEband(n, A1d, b, temp, nrhs)
  DO j = 1,n
    DO i = 1,n
      b(j,i) = temp(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE matvec2d


SUBROUTINE shftmat2dumb(shftmat,hmat,omat,a,nb1,nb2,ldab)
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
    WRITE(*,*) 'fatal error in shftmat2dumb: ldab (=', ldab, ') is too small'
    STOP 181
  ENDIF
  centroid = kl + ku + 1
  strd1 = 2*d1 + 1
  strd2 = MAX(nb1,strd1)
  shftmat(:,1:nb1*nb2) = 0.d0
  DO col = 1, nb1*nb2
    DO blk = -d2,d2
      id = blk*strd2 + centroid
      jd = blk*strd1
      shftmat(id-d1:id+d1,col) = hmat(jd-d1:jd+d1,col) - a*omat(jd-d1:jd+d1,col)
    END DO
  END DO

  RETURN
END SUBROUTINE shftmat2dumb


SUBROUTINE dumpbandmatrix(A,ldab,n,kl,ku)
  implicit none

  INTEGER :: ldab, n, kl, ku
  REAL(kind=8) :: A(ldab,n)
  INTEGER :: row, col, j, shft
  REAL(kind=8) :: tol

  tol = 1.d-50
  IF ( ldab >= 2*kl + ku + 1 ) THEN
    shft = kl + ku + 1
  ELSEIF ( ldab >= kl + ku + 1 ) THEN
    shft = ku + 1
  ELSE
    WRITE(*,*) 'ERROR in dumpbandmatrix: kl, ku, and ldab are incompatible'
    STOP 151
  ENDIF
  OPEN(UNIT=8,FILE='A.dat')
  DO col = 1,n
    DO j = 1,ldab
      IF (abs(A(j,col)) > tol) THEN
        row = j + col - shft
        WRITE(8,100) row, col, A(j,col)
      ENDIF
    END DO
  END DO
  CLOSE(8)
  WRITE(*,*) 'wrote matrix to A.dat'

  RETURN
100 FORMAT(2I4,E14.6)
END SUBROUTINE dumpbandmatrix
