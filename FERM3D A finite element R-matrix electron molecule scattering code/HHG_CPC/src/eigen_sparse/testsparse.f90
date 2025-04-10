!!!*************************************************************
! 文件/File: testsparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: testsparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM testsparse
  USE FiniteElements
  USE potentials
  implicit none

  REAL(kind=8), EXTERNAL :: efield

  TYPE(FEBnodes) :: x
  INTEGER, ALLOCATABLE :: rowind(:), colptr(:)
  REAL(kind=8), ALLOCATABLE :: FEmat(:,:), Hmat(:), Omat(:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:), cmask(:)
  REAL(kind=8), ALLOCATABLE :: psico(:,:), workK(:,:)
  CHARACTER ::  pot_type*10
  INTEGER :: nmnds, nbasis, nnz, nsave, nprop
  INTEGER :: nev, j, nconv, ldab, info, coarsesteps, nrhs
  INTEGER, ALLOCATABLE :: ipiv(:)
  REAL(kind=8) :: xmax, sigmar, t0, E0, savetime, t, dt, norm, imag, Amp, alpha, beta

  t0 = 0.d0
  E0 = 1.d-3
  nsave = 570
  savetime = 1.d0
  nprop = 100

  CALL CreateNodes(x,'nodes.inp')
  !!CALL PrintNodes(x,6)
  nmnds = x%numnds

  nbasis = 3*(nmnds-2) + 4
  IF (nmnds <= 2) THEN
    nnz = nbasis*nbasis
  ELSE
    nnz = x%numleft*(x%numleft+6) + x%numright*(x%numright+6) + &
        & 9*(3*nmnds-8)
  ENDIF
  ALLOCATE( FEmat(-5:5,nbasis), Hmat(nnz), Omat(nnz), rowind(nnz), colptr(0:nbasis) )
  pot_type = 'coulomb'
  alpha = -1.d0
  beta = 2.d0
  CALL Vmtrx(FEmat,pot_type,alpha,beta,x)
  CALL sparse1d(Omat,rowind,colptr,nbasis,nnz,FEmat)
  CALL KEmtrx(FEmat,x)
  CALL sparsegrab1d(Hmat,rowind,colptr,nbasis,nnz,FEmat)
  Hmat(:) = Hmat(:) + Omat(:)
  CALL Omtrx(FEmat,x)
  CALL sparsegrab1d(Omat,rowind,colptr,nbasis,nnz,FEmat)

  sigmar = -0.7d0
  nev   = 6 
  ALLOCATE( eval(nev), evec(nbasis,nev) )
  eval(:) = 0.d0
  evec(:,:) = 0.d0
  CALL eigensparse(Hmat,Omat,rowind,colptr,nbasis,nnz,sigmar,eval,evec,nev,nconv,4)

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues--'
  DO j = 1,nev
    WRITE(*,*) j, ':  ', eval(j)
  END DO

  WRITE(*,*)
  STOP
END PROGRAM testsparse
