!!!*************************************************************
! 文件/File: prop2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: prop2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM prop2d
  USE FiniteElements
  implicit none

  REAL(kind=8), EXTERNAL :: efield

  TYPE(FEBnodes) :: x
  REAL(kind=8), ALLOCATABLE :: H0mat(:), Omat(:), Hmat(:), EFmat(:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:)
  REAL(kind=8), ALLOCATABLE :: psico(:,:), workK(:,:), cmask(:)
  INTEGER :: nsave, nprop, nnz, nb2
  INTEGER :: nev, j, nconv, coarsesteps, col, nrhs, info
  INTEGER, ALLOCATABLE :: rowind(:), colptr(:), rowind_save(:)
  REAL(kind=8) :: xmax, sigmar, t0, E0, savetime, t, dt, norm, imag, Amp
  CHARACTER :: save_fact

  t0 = 0.d0
  E0 = 1.d-3
  nsave = 10
  savetime = 1.d0
  nprop = 100

  CALL CreateNodes(x,'nodes.inp')
  CALL PrintNodes(x,6)

  nb2 = (x%numbasis)**2
  IF (x%numnds <= 2) THEN
    nnz = nb2
  ELSE
    nnz = x%numleft*(x%numleft+6) + x%numright*(x%numright+6) + &
	& 9*(3*x%numnds-8)
  ENDIF
  nnz = nnz*nnz
  WRITE(*,*) 'guessed number of nonzeros is ', nnz
  ALLOCATE( H0mat(nnz), rowind(nnz), colptr(0:nb2), Omat(nnz) )
  CALL creatematrices(H0mat, Omat, rowind, colptr, nnz, x)
  WRITE(*,*) 'actual number of nonzeros is ', nnz

  sigmar = -1.0d0
  nev   = 6 
  ALLOCATE( eval(nev), evec(nb2,nev) )
  eval = 0.d0
  evec = 0.d0
  CALL eigensparse(H0mat,Omat,rowind,colptr,nb2,nnz,sigmar,eval,evec,nev,nconv,0)

  ALLOCATE( psico(nb2,2), cmask(nb2), workK(nb2,10) )
  workK(:,:) = 0.d0
  psico(:,1) = evec(:,1)
  psico(:,2) = 0.d0
  CALL coefmask(cmask,x,4,4)
  WRITE(*,*) 'initial state energy = ', eval(1)
  CALL sparsemultiply(workK(1,1),psico,Omat,rowind,colptr,nb2,nnz,2)
  norm = DOT_PRODUCT(psico(:,1),workK(:,1)) + DOT_PRODUCT(psico(:,2),workK(:,2))
  WRITE(*,*) 'norm of initial state is ', norm

  DEALLOCATE( eval, evec )

  ALLOCATE( rowind_save(nnz), EFmat(nnz) )
  !! save original overlap matrix
  EFmat(:) = Omat(:)
  rowind_save(:) = rowind(:)
  ! workK used as dummy variable, not used
  save_fact = 'Y'
  nrhs = 0
  ! LU factorize the overlap matrix
  CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK,1,info,save_fact,2)
  rowind(:) = rowind_save(:)
  Omat(:) = EFmat(:)
  DEALLOCATE( rowind_save )

  dt = savetime/nprop
  t = t0
  nrhs = 2 !!! real/complex parts will be stored in separate columns
  ALLOCATE( Hmat(nnz) )
  Hmat(:) = H0mat(:)

  DO coarsesteps = 1,nsave
  DO j = 1,nprop

    !TTT!Amp = E0*efield(t-t0)
    !TTT!Hmat = H0mat + Amp*EFmat
!! build K1 into workK(:,1:2)
    CALL sparsemultiply(workK(1,1),psico,Hmat,rowind,colptr,nb2,nnz,nrhs)
    CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK(1,1), &
                 & nb2,info,save_fact,2)
    !! Re(K1) = workK(:,2),  Im(K1) = workK(:,1)
    workK(:,1) = -dt*workK(:,1)
    workK(:,2) = dt*workk(:,2)
   
    !TTT!Amp = E0*efield(t+dt/2-t0)
    !TTT!Hmat = H0mat + Amp*EFmat
!! build K2 into workK(:,3:4), use workK(:,5:6) as scratch
    workK(:,5) = psico(:,1) + 0.5d0*workK(:,2)
    workK(:,6) = psico(:,2) + 0.5d0*workK(:,1)
    CALL sparsemultiply(workK(1,3),workK(1,5),Hmat,rowind,colptr,nb2,nnz,nrhs)
    CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK(1,3),nb2, &
                 & info,save_fact,2)
    !! Re(K2) = workK(:,4),  Im(K2) = workK(:,3)
    workK(:,3) = -dt*workK(:,3)
    workK(:,4) = dt*workK(:,4)

!! build K3 into workK(:,5:6), use workK(:,7:8) as scratch
    workK(:,7) = psico(:,1) + 0.5d0*workK(:,4)
    workK(:,8) = psico(:,2) + 0.5d0*workK(:,3)
    CALL sparsemultiply(workK(1,5),workK(1,7),Hmat,rowind,colptr,nb2,nnz,nrhs)
    CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK(1,5),nb2, &
                 & info,save_fact,2)
    !! Re(K3) = workK(:,6),  Im(K3) = workK(:,5)
    workK(:,5) = -dt*workK(:,5)
    workK(:,6) = dt*workK(:,6)

    !TTT!Amp = E0*efield(t+dt-t0)
    !TTT!Hmat = H0mat + Amp*EFmat
!! build K4 into workK(:,7:8), use workK(:,9:10) as scratch
    workK(:,9) = psico(:,1) + workK(:,6)
    workK(:,10) = psico(:,2) + workK(:,5)
    CALL sparsemultiply(workK(1,7),workK(1,9),Hmat,rowind,colptr,nb2,nnz,nrhs)
    CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK(1,7),nb2, &
                 & info,save_fact,2)
    !! Re(K4) = workK(:,8),  Im(K4) = workK(:,7)
    workK(:,7) = -dt*workK(:,7)
    workK(:,8) = dt*workK(:,8)

    psico(:,1) = psico(:,1) + (workK(:,2)+2*workK(:,4)+2*workK(:,6)+workK(:,8))/6.d0
    psico(:,2) = psico(:,2) + (workK(:,1)+2*workK(:,3)+2*workK(:,5)+workK(:,7))/6.d0
    !TTT!psico(:,1) = cmask(:)*psico(:,1)
    !TTT!psico(:,2) = cmask(:)*psico(:,2)
    t = t + dt

  END DO

  nrhs = 2
  CALL sparsemultiply(workK(1,9),psico,Omat,rowind,colptr,nb2,nnz,nrhs)
  norm = DOT_PRODUCT(psico(:,1),workK(:,9)) + DOT_PRODUCT(psico(:,2),workK(:,10))
  imag = DOT_PRODUCT(psico(:,1),workK(:,10)) - DOT_PRODUCT(psico(:,2),workK(:,9))
  WRITE(*,*) t, ':  dnorm = ', 1.d0 - norm
  END DO

  !! clear LU factorization
  save_fact = 'N'
  nrhs = 0
  CALL dSLUsolve(nb2,nnz,nrhs,Omat,rowind,colptr,workK,1,info,save_fact,0)
  DEALLOCATE( EFmat )

  STOP
END PROGRAM prop2d



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
