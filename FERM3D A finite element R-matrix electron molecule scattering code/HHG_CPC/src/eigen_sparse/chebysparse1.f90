!!!*************************************************************
! 文件/File: chebysparse1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: chebysparse1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


! psico(nb2,2): real part stored in first column, imaginary in second column
!               on entry, initial state coefficients
!               on exit, final state coefficients
! nb2: dimension of psico and matrices
! spect(2): bounds of eigenvalues
!           spect(1) = lower bound, spect(2) = upper bound

SUBROUTINE chebysparse1(psico,dt,Hmat,rowind,colptr,nb2,nnz,spect)
  implicit none

  INTEGER :: nb2, nnz
  INTEGER :: rowind(nnz), colptr(0:nb2)
  REAL(kind=8) :: Hmat(nnz), psico(nb2,2)
  REAL(kind=8) :: dt, spect(2)
  INTEGER :: maxnmax, nmax, nrhs, info, n, j
  REAL(kind=8), ALLOCATABLE :: bessj(:), Wmat(:), Tn(:,:), Tnm1(:,:), Tnm2(:,:)
  REAL(kind=8) :: fco, gco, phase(2)


  fco = 0.5d0*dt*(spect(2)-spect(1))
  gco = spect(1)*dt
  phase(1) = cos( fco + gco )
  phase(2) = -sin( fco + gco )
  !WRITE(*,*) 'fco = ', fco

  maxnmax = (3.4d0/abs(fco) + 4.2d0/sqrt(abs(fco)) + 1.00d0)*abs(fco) + 2
  maxnmax = max(maxnmax,7)
  !WRITE(*,*) 'maxnmax = ', maxnmax
  ALLOCATE( bessj(0:maxnmax) )
  bessj(:) = 0.d0
  CALL makebessel(bessj, fco, maxnmax)
  nmax = 2
  DO n = maxnmax,1,-1
   IF ( abs(bessj(n)) > 3.d-10 ) THEN
     nmax = n
     EXIT
   ENDIF
  END DO
  nmax = MIN(MAX(nmax,5),maxnmax)
  !WRITE(*,*) '   nmax = ', nmax
  !WRITE(*,*) 'last coefficient = ', bessj(nmax)
  !WRITE(*,*)

  ALLOCATE( Wmat(nnz) )

  !WRITE(*,*) 'dt = ', dt
  !WRITE(*,*) 'fco = ', fco
  !WRITE(*,*) 'gco = ', gco
  Wmat(:) = (-dt/fco)*Hmat(:) ! + (1.d0 + gco/fco)*Omat(:)
  DO n = 0,(nb2-1)
   DO j = colptr(n)+1,colptr(n+1)
     IF ( rowind(j) == n) THEN
       Wmat(j) = Wmat(j) + (1.d0 + gco/fco)
       !WRITE(*,*) rowind(j)+1, n+1, nb2
     ENDIF
   END DO
  END DO

  ALLOCATE( Tn(nb2,2), Tnm1(nb2,2), Tnm2(nb2,2) )
  Tn(:,:) = 0.d0
  Tnm1(:,:) = 0.d0
  Tnm2(:,:) = 0.d0
  nrhs = 2 !!! real/complex parts will be stored in separate columns

  IF (nmax >= 0) THEN
!! n = 0 term
    Tn(:,:) = psico(:,:)
    psico(:,1) = bessj(0)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
    psico(:,2) = bessj(0)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
  ENDIF

!! n = 1 term
  IF (nmax >= 1) THEN
    Tnm1(:,:) = Tn(:,:)
    CALL sparsemultiply(Tn,Tnm1,Wmat,rowind,colptr,nb2,nnz,nrhs)
    psico(:,1) = psico(:,1) - 2*bessj(1)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
    psico(:,2) = psico(:,2) + 2*bessj(1)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
  ENDIF

!! n > 1 terms
  DO n = 2,nmax
    Tnm2(:,:) = Tnm1(:,:)
    Tnm1(:,:) = Tn(:,:)
    CALL sparsemultiply(Tn,Tnm1,Wmat,rowind,colptr,nb2,nnz,nrhs)
    Tn(:,:) = 2*Tn(:,:) - Tnm2(:,:)
    SELECT CASE ( MOD(n,4) )
    CASE(0)
      psico(:,1) = psico(:,1) + 2*bessj(n)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
      psico(:,2) = psico(:,2) + 2*bessj(n)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
    CASE(1)
      psico(:,1) = psico(:,1) - 2*bessj(n)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
      psico(:,2) = psico(:,2) + 2*bessj(n)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
    CASE(2)
      psico(:,1) = psico(:,1) - 2*bessj(n)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
      psico(:,2) = psico(:,2) - 2*bessj(n)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
    CASE(3)
      psico(:,1) = psico(:,1) + 2*bessj(n)*(phase(1)*Tn(:,2) + phase(2)*Tn(:,1))
      psico(:,2) = psico(:,2) - 2*bessj(n)*(phase(1)*Tn(:,1) - phase(2)*Tn(:,2))
    END SELECT
  END DO

  RETURN
END SUBROUTINE chebysparse1

