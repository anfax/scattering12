!!!*************************************************************
! 文件/File: fd2.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: fd2.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


PROGRAM fd2
  implicit none

  INTEGER :: j, l, count, nsave, nprop
  INTEGER :: ngrd, nmat, nnz, mid, nev, nconv, nlftmsk, nrgtmsk
  INTEGER, ALLOCATABLE :: rowind(:), colptr(:)
  REAL(kind=8) :: dx, xmax, shift, spect0(2), spect(2), dt, deltat, E0, amp, savetime, t
  REAL(kind=8), ALLOCATABLE :: x(:), V1d(:), Vefield(:), lftmsk(:), rgtmsk(:)
  REAL(kind=8), ALLOCATABLE :: psi(:,:), psireal(:,:), psiimag(:,:)
  REAL(kind=8), ALLOCATABLE :: Hmat(:), H0mat(:), H1d(:,:), eval(:), evec(:,:)

  dx = 1.d-1
  xmax = 30
  ngrd = 2*INT(xmax/dx) + 1
  ALLOCATE( x(ngrd), V1d(ngrd), H1d(-3:3,ngrd) )
  x(:) = 0.d0
  V1d(:) = 0.d0
  H1d(:,:) = 0.d0
  mid = ngrd/2 + MOD(ngrd,2)
  DO j = 1, (ngrd-mid)
    x(mid+j) = j*dx
    x(mid-j) = -j*dx
  END DO
  xmax = x(ngrd)
  V1d(:) = -1.d0/sqrt(x(:)*x(:) + 2.d0)

  H1d(-3,1) = -0.5d0/(90.d0*dx*dx)
  H1d(-2,1) = -13.5d0*H1d(-3,1)
  H1d(-1,1) = -10.d0*H1d(-2,1)
  H1d(0,1) = -245.d0*H1d(-3,1)
  H1d(1,1) = H1d(-1,1)
  H1d(2,1) = H1d(-2,1)
  H1d(3,1) = H1d(-3,1)
  DO j = 2,ngrd
    H1d(:,j) = H1d(:,1)
  END DO
  DO j = 1,3
    H1d(-3:-j,j) = 0.d0
    H1d(j:3,ngrd+1-j) = 0.d0
  END DO
  H1d(0,:) = H1d(0,:) + V1d(:)

  nnz = 7*ngrd - 6 - 6
  nmat = ngrd
  !nnz = nnz*nnz
  !nmat = nmat*nmat
  ALLOCATE( rowind(nnz), colptr(0:nmat), Hmat(nnz), H0mat(nnz) )
  rowind(:) = 0
  colptr(:) = 0
  count = 0
  DO j = 1,3
    DO l = (1-j),3
      count = count + 1
      colptr(j) = colptr(j) + 1
      rowind(count) = j + l - 1
      H0mat(count) = H1d(-l,j+l)
    END DO
  END DO
  DO j = 4,(ngrd-3)
    DO l = -3,3
      count = count + 1
      colptr(j) = colptr(j) + 1
      rowind(count) = j + l - 1
      H0mat(count) = H1d(-l,j+l)
    END DO
  END DO
  DO j = (ngrd-2),ngrd
    DO l = -3,ngrd-j
      count = count + 1
      colptr(j) = colptr(j) + 1
      rowind(count) = j + l - 1
      H0mat(count) = H1d(-l,j+l)
    END DO
  END DO

  DO j = 2,nmat
    colptr(j) = colptr(j) + colptr(j-1)
  END DO
  WRITE(*,*) 'check:  ', count, nnz

!---- add Efield potential
  E0 = 1.0d-1
  spect0(1) = -1.d0
  spect0(2) = 3.142d0/dx/dx
  spect(1) = spect0(1) - E0*xmax
  spect(2) = spect0(2) + E0*xmax
  WRITE(*,*) 'estimated max eigenvalue is ', spect(2)

  ALLOCATE( Vefield(nmat) )
  Vefield(:) = E0*x(:)
  Hmat(:) = H0mat(:)

  shift = -1.0d0
  nev = 3
  ALLOCATE( eval(nev), evec(nmat,nev) )
  CALL egnsps1(H0mat,rowind,colptr,nmat,nnz,shift,eval,evec,nev,nconv,2)

  WRITE(*,*)
  WRITE(*,*) 'computed ', nconv, ' eigenvalues'
  DO j = 1,nconv
    WRITE(*,*) j, eval(j)
  END DO

  nlftmsk = 20.0d0/dx
  nrgtmsk = nlftmsk
  ALLOCATE( lftmsk(nlftmsk), rgtmsk(nrgtmsk) )
  DO j = 1,nlftmsk
    lftmsk(j) = mask(x(j),x(nlftmsk),x(1))
  END DO
  DO j = 1,nlftmsk
    rgtmsk(j) = lftmsk(nlftmsk+1-j)
  END DO

  deltat = 1.d0
  nsave = 6.283d0*50/0.55d0/deltat

  ALLOCATE( psi(nmat,2), psireal(nmat,nsave), psiimag(nmat,nsave) )
  psi(:,1) = evec(:,1)
  psi(:,2) = 0.d0
  nprop = 20
  dt = deltat/nprop
  savetime = 0.d0

  DO j = 1,nsave
    t = savetime + 0.5*dt
    DO count = 1,nprop
      amp = efield(t)
      CALL addefield(amp)
      CALL chebysparse1(psi,dt,Hmat,rowind,colptr,nmat,nnz,spect)
      psi(1:nlftmsk,1) = lftmsk(:)*psi(1:nlftmsk,1)
      psi(1:nlftmsk,2) = lftmsk(:)*psi(1:nlftmsk,2)
      psi(nmat-nrgtmsk+1:nmat,1) = rgtmsk(:)*psi(nmat-nrgtmsk+1:nmat,1)
      psi(nmat-nrgtmsk+1:nmat,2) = rgtmsk(:)*psi(nmat-nrgtmsk+1:nmat,2)
      t = t + dt
    END DO
    psireal(:,j) = psi(:,1)
    psiimag(:,j) = psi(:,2)
    savetime = savetime + deltat
    WRITE(*,*) savetime, amp
  END DO

  OPEN(UNIT=8, FILE='fd2.bin', FORM='UNFORMATTED')
  WRITE(8) ngrd, x, nmat, nsave, psireal, psiimag
  CLOSE(8)
  WRITE(*,*) 'saved data to fd2.bin'


100 FORMAT(I7,I7,E13.5)
101 FORMAT(2E13.5)

  STOP

  CONTAINS

  REAL(kind=8) FUNCTION mask(xx,ndi,ndx)
    REAL(kind=8) :: xx, ndi, ndx, zx, prec, ac, dz

    prec = 0.01d0
    dz = abs(ndx-ndi)
    IF ( dz < 1.d-40 ) THEN
      mask = 1.d0
    ELSE
      zx = abs(ndx-xx)
      ac = -2.d0*log(2.d0*prec)/log(2.d0-2.d0*prec)/dz
      mask = 0.5d0*(1.d0+ tanh(ac*(zx-0.5d0*dz)))
    ENDIF

    !prec = 0.0001d0
    !dz = abs(ndx-ndi)
    !zx = abs(ndx-xx)
    !ac = -2.d0*log(2.d0*prec)/log(2.d0-2.d0*prec)/dz
    !mask = 0.5d0*(1.d0+ tanh(ac*(zx-0.5d0*dz)))

    RETURN
  END FUNCTION mask


  REAL(kind=8) FUNCTION efield(tau)
    REAL(kind=8) :: tau
    REAL(kind=8) :: omega, cycles
  
    omega = 0.55d0
    cycles = 50.d0
    efield = ((sin(0.5*omega*tau/cycles))**2)*sin(omega*tau)

    RETURN
  END FUNCTION efield


  SUBROUTINE addefield(amp)
    REAL(kind=8) :: amp
    INTEGER :: l, j

    DO l = 0,(nmat-1)
      DO j = colptr(l)+1,colptr(l+1)
        IF ( rowind(j) == l) THEN
          Hmat(j) = H0mat(j) + amp*Vefield(l+1)
        ENDIF
      END DO
    END DO

    RETURN
  END SUBROUTINE addefield

END PROGRAM fd2
