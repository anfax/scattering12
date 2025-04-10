!!!*************************************************************
! 文件/File: cheby2.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: cheby2.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM cheby2
  USE FiniteElements
  implicit none

  REAL(kind=8), EXTERNAL :: efield

  TYPE(FEBnodes) :: x
  REAL(kind=8), ALLOCATABLE :: Hmat(:), Omat(:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:)
  REAL(kind=8), ALLOCATABLE :: psico(:,:), psi0(:,:)
  INTEGER :: nnz, nb2
  INTEGER :: nev, n, j, nconv, col, nrhs, info
  INTEGER, ALLOCATABLE :: rowind(:), colptr(:)
  REAL(kind=8) :: xmax, sigmar, norm, lamdamax, lamdamin, pi, deltat, spect(2)
  CHARACTER :: save_fact

  WRITE(*,*) 'enter time step'
  READ(*,*) deltat
  pi = 4.d0*atan(1.d0)

  CALL CreateNodes(x,'nodes.inp')
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

  nb2 = (x%numbasis)**2
  IF (x%numnds <= 2) THEN
    nnz = nb2
  ELSE
    nnz = x%numleft*(x%numleft+6) + x%numright*(x%numright+6) + &
	& 9*(3*x%numnds-8)
  ENDIF
  nnz = nnz*nnz
  WRITE(*,*) 'guessed number of nonzeros is ', nnz
  ALLOCATE( Hmat(nnz), rowind(nnz), colptr(0:nb2), Omat(nnz) )
  CALL creatematrices(Hmat, Omat, rowind, colptr, nnz, x)
  WRITE(*,*) 'actual number of nonzeros is ', nnz

  sigmar = -1.0d0
  nev   = 6 
  ALLOCATE( eval(nev), evec(nb2,nev) )
  eval = 0.d0
  evec = 0.d0
  CALL eigensparse(Hmat,Omat,rowind,colptr,nb2,nnz,sigmar,eval,evec,nev,nconv,0)

  ALLOCATE( psico(nb2,2), psi0(nb2,2) )
  psi0(:,:) = 0.d0
  psico(:,1) = evec(:,1)
  psico(:,2) = 0.d0
  WRITE(*,*) 'initial state energy = ', eval(1)
  CALL sparsemultiply(psi0(1,1),psico,Omat,rowind,colptr,nb2,nnz,2)
  norm = DOT_PRODUCT(psico(:,1),psi0(:,1)) + DOT_PRODUCT(psico(:,2),psi0(:,2))
  WRITE(*,*) 'norm of initial state is ', norm
  !! save initial state
  psi0(:,:) = psico(:,:)

  DEALLOCATE( evec )

  CALL chebysparse(psico,deltat,Hmat,Omat,rowind,colptr,nb2,nnz,spect,'N')


  nrhs = 2
  !CALL sparsemultiply(psi0(1,1),psico,Omat,rowind,colptr,nb2,nnz,nrhs)
  !norm = DOT_PRODUCT(psico(:,1),psi0(:,1)) + DOT_PRODUCT(psico(:,2),psi0(:,2))
  !WRITE(*,*) ' norm = ', norm

  psi0(:,2) = -sin(deltat*eval(1))*psi0(:,1) - psico(:,2)
  psi0(:,1) = cos(deltat*eval(1))*psi0(:,1) - psico(:,1)
  norm = DOT_PRODUCT(psi0(:,1),psi0(:,1)) + DOT_PRODUCT(psi0(:,2),psi0(:,2))
  WRITE(*,*) 'error = ', norm

  STOP
END PROGRAM cheby2



SUBROUTINE creatematrices(Hsparse, Osparse, rowind, colptr, nnzinout, x)
  USE FiniteElements
  USE Potentials
  implicit none

  INTEGER :: nnzinout
  INTEGER :: rowind(nnzinout), colptr(0:*)
  REAL(kind=8) :: Hsparse(nnzinout), Osparse(nnzinout)
  TYPE(FEBnodes) :: x
  INTEGER :: nbasis, nnz, nb2
  REAL(kind=8), ALLOCATABLE :: O1d(:,:), H1d(:,:), V1d(:,:), H2d(:,:), O2d(:,:)
  REAL(kind=8) :: alpha, beta
  CHARACTER :: pot1_type*10, pot2_type*10

  pot1_type = 'coulomb'
  alpha = -1.d0
  beta = 2.d0
  pot2_type = 'zero'

  nbasis = x%numbasis
  nb2 = nbasis*nbasis
  ALLOCATE( O1d(-5:5,nbasis), H1d(-5:5,nbasis), V1d(-5:5,nbasis) )
  CALL Vmtrx(V1d,pot1_type,alpha,beta,x)
  CALL KEmtrx(H1d,x)
  H1d(:,:) = H1d(:,:) + V1d(:,:)
  CALL Omtrx(O1d,x)

  ALLOCATE( H2d(-60:60,nb2), O2d(-60:60,nb2) )
  H2d(:,:) = 0.d0
  O2d(:,:) = 0.d0

  CALL ham2d(H2d,H1d,O1d,H1d,O1d,nbasis,nbasis)
  !! O2d temporarily holds the 2-electron potential
  CALL V2dmtrx(O2d,pot2_type,x,x)
  H2d(:,:) = H2d(:,:) + O2d(:,:)

  CALL directsum(O2d,O1d,O1d,nbasis,nbasis)

  CALL sparse2d(Hsparse,rowind,colptr,nb2,nnz,H2d,nbasis,nbasis)

  CALL sparsegrab2d(Osparse,rowind,colptr,nb2,nnz,O2d,nbasis,nbasis)

  nnzinout = nnz
  RETURN
END SUBROUTINE creatematrices

