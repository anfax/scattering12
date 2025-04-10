!!!*************************************************************
! 文件/File: prop1d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: prop1d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

PROGRAM prop1d
  USE FiniteElements
  USE potentials
  implicit none

  REAL(kind=8), EXTERNAL :: efield

  TYPE(FEBnodes) :: x
  REAL(kind=8), ALLOCATABLE :: Omat(:,:), Ofact(:,:), H0mat(:,:), Vmat(:,:), Hmat(:,:)
  REAL(kind=8), ALLOCATABLE :: eval(:), evec(:,:), cmask(:)
  REAL(kind=8), ALLOCATABLE :: psico(:,:), workK(:,:)
  CHARACTER ::  pot_type*10
  INTEGER :: nmnds, nbasis, nsave, nprop
  INTEGER :: nev, j, nconv, ldab, info, ierr, coarsesteps, nrhs
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
  ALLOCATE( Omat(-5:5,nbasis), H0mat(-5:5,nbasis), Vmat(-5:5,nbasis) )
  pot_type = 'coulomb'
  alpha = -1.d0
  beta = 2.d0
  CALL Vmtrx(Vmat,pot_type,alpha,beta,x)
  CALL KEmtrx(H0mat,x)
  H0mat(:,:) = H0mat(:,:) + Vmat(:,:)
  CALL Omtrx(Omat,x)

  ALLOCATE( eval(nbasis), cmask(nbasis) )
  eval = 0.d0
  cmask = 1.d0
  CALL Av_FEband(nbasis, Omat, cmask, eval, 1)
  DEALLOCATE( eval, cmask )

  sigmar = -0.7d0
  nev   = 6 
  ALLOCATE( eval(nev), evec(nbasis,nev) )
  CALL eigenfew(H0mat,Omat,nbasis,sigmar,eval,evec,nev,nconv,0)


  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues--'
  DO j = 1,nev
    WRITE(*,*) j, ':  ', eval(j)
  END DO

  STOP

  pot_type = 'power'
  alpha = 1.d0
  beta = 1.d0
  CALL Vmtrx(Vmat,pot_type,alpha,beta,x)

  ldab = 2*5 + 5 + 1
  ALLOCATE( Ofact(ldab,nbasis), ipiv(nbasis) )
  DO j = 1,nbasis
    Ofact(6:ldab,j) = Omat(-5:5,j)
  END DO
  CALL dgbtrf(nbasis,nbasis,5,5,Ofact,ldab,ipiv,info)

!! psico and workK are complex, with real/imaginary in adjacent columns
  ALLOCATE( Hmat(-5:5,nbasis) )
  ALLOCATE( psico(nbasis,2), cmask(nbasis), workK(nbasis,10) )
  psico(:,1) = evec(:,1)
  psico(:,2) = 0.d0
  CALL coefmask(cmask,x,4,4)

  dt = savetime/nprop
  t = t0
  nrhs = 2 !!! real/complex parts will be stored in separate columns

  DO coarsesteps = 1,nsave
  DO j = 1,nprop

    Amp = E0*efield(t-t0)
    Hmat = H0mat + Amp*Vmat
!! build K1 into workK(:,1:2)
    CALL Av_FEband(nbasis, Hmat, psico, workK(1,1), nrhs)
    CALL dgbtrs('N',nbasis,5,5,nrhs,Ofact,ldab,ipiv,workK(1,1),nbasis,ierr)
    !! Re(K1) = workK(:,2),  Im(K1) = workK(:,1)
    workK(:,1) = -dt*workK(:,1)
    workK(:,2) = dt*workk(:,2)
    
    Amp = E0*efield(t+dt/2-t0)
    Hmat = H0mat + Amp*Vmat
!! build K2 into workK(:,3:4), use workK(:,5:6) as scratch
    workK(:,5) = psico(:,1) + 0.5d0*workK(:,2)
    workK(:,6) = psico(:,2) + 0.5d0*workK(:,1)
    CALL Av_FEband(nbasis, Hmat, workK(1,5), workK(1,3), nrhs)
    CALL dgbtrs('N',nbasis,5,5,nrhs,Ofact,ldab,ipiv,workK(1,3),nbasis,ierr)
    !! Re(K2) = workK(:,4),  Im(K2) = workK(:,3)
    workK(:,3) = -dt*workK(:,3)
    workK(:,4) = dt*workK(:,4)
    
!! build K3 into workK(:,5:6), use workK(:,7:8) as scratch
    workK(:,7) = psico(:,1) + 0.5d0*workK(:,4)
    workK(:,8) = psico(:,2) + 0.5d0*workK(:,3)
    CALL Av_FEband(nbasis, Hmat, workK(1,7), workK(1,5), nrhs)
    CALL dgbtrs('N',nbasis,5,5,nrhs,Ofact,ldab,ipiv,workK(1,5),nbasis,ierr)
    !! Re(K3) = workK(:,6),  Im(K3) = workK(:,5)
    workK(:,5) = -dt*workK(:,5)
    workK(:,6) = dt*workK(:,6)

    Amp = E0*efield(t+dt-t0)
    Hmat = H0mat + Amp*Vmat
!! build K4 into workK(:,7:8), use workK(:,9:10) as scratch
    workK(:,9) = psico(:,1) + workK(:,6)
    workK(:,10) = psico(:,2) + workK(:,5)
    CALL Av_FEband(nbasis, Hmat, workK(1,9), workK(1,7), nrhs)
    CALL dgbtrs('N',nbasis,5,5,nrhs,Ofact,ldab,ipiv,workK(1,7),nbasis,ierr)
    !! Re(K4) = workK(:,8),  Im(K4) = workK(:,7)
    workK(:,7) = -dt*workK(:,7)
    workK(:,8) = dt*workK(:,8)

    psico(:,1) = psico(:,1) + (workK(:,2)+2*workK(:,4)+2*workK(:,6)+workK(:,8))/6.d0
    psico(:,2) = psico(:,2) + (workK(:,1)+2*workK(:,3)+2*workK(:,5)+workK(:,7))/6.d0
    psico(:,1) = cmask(:)*psico(:,1)
    psico(:,2) = cmask(:)*psico(:,2)
    t = t + dt

  END DO
  CALL Av_FEband(nbasis, Omat, psico, workK(1,9), nrhs)
  norm = DOT_PRODUCT(psico(:,1),workK(:,9)) + DOT_PRODUCT(psico(:,2),workK(:,10))
  imag = DOT_PRODUCT(psico(:,1),workK(:,10)) - DOT_PRODUCT(psico(:,2),workK(:,9))
  WRITE(*,*) t, ':  dnorm = ', 1.d0 - norm
  END DO

  STOP
END PROGRAM prop1d
 
