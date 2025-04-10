!!!*************************************************************
! 文件/File: lin_3drmatsparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: lin_3drmatsparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

module Solve
  use nrtype, only : dbl,i4b
  real(kind=dbl),dimension(:),allocatable :: eigen_value
  real(kind=dbl),dimension(:,:),allocatable :: eigen_vector 
  real(kind=dbl) :: Z
  integer(kind=i4b) :: rowA

  ! interface invert
  !   module procedure invert_real,invert_complex
  ! end interface invert

  !----------------------------
  ! This module contains all the subroutines that perform the linear algebra
  ! tasks or setup them
  !----------------------------

contains


  subroutine lin_solve(a,b,n,lda,nrhs)
    use nrtype, only : dbl,i4b
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,nrhs,n,m
    character trans
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(n,1) :: B
    integer(kind=i4b),dimension(:),allocatable :: ipiv 
    !------------------------------------------------------------------------
    ! Subroutine used to solve the linear system that rises in Laplace/Poisson equation
    ! Electrostatic potential problems
    !------------------------------------------------------------------------

    lda=n
    trans='n'
    m=lda
    ldb=lda 
    allocate(ipiv(n)) 
    call dgetrf(m,n,a,lda,ipiv,info)
    call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
90  format(4d14.7)
    deallocate(ipiv)
  end subroutine lin_solve


  subroutine invert_real(a,n,lda)
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,n,m,lwork
    real(kind=dbl),dimension(n,n) :: A
    integer(kind=i4b),dimension(:),allocatable :: ipiv 
    real(kind=dbl),dimension(:),allocatable :: work

    !------------------------------------------------------------------------
    ! Inverts the matrix omega
    !------------------------------------------------------------------------
    ! The dimensioning of lwork is not the best possible
    !------------------------------------------------------------------------
    lda=n
    m=lda
    ldb=lda 
    lwork=4*n
    allocate(ipiv(n)) 
    allocate(work(lwork)) 
    write(6,*)m,n,lda
    call dgetrf(m,n,a,lda,ipiv,info)
!    write(6,*)'after dgetrf'
    call dgetri(n,a,lda,ipiv,work,lwork,info)
!    write(6,*)'after dgetrs  info=',info
90  format(4d14.7)
    deallocate(ipiv)
    deallocate(work)
    write(6,*)'after invert real'
  end subroutine invert_real



  subroutine invert_complex(a,n,lda)
    use nrtype,only : dbl,i4b,dpc
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,n,m,lwork
    complex(kind=dpc),dimension(n,n) :: A
    integer(kind=i4b),dimension(:),allocatable :: ipiv
    complex(kind=dpc),dimension(:),allocatable :: work
    !------------------------------------------------------------------------
    ! Inverts the matrix omega
    !------------------------------------------------------------------------
    ! The dimensioning of lwork is not the best possible
    !------------------------------------------------------------------------
    lda=n
    m=lda
    ldb=lda
    lwork=4*n
    allocate(ipiv(n))
    allocate(work(lwork))
    call zgetrf(m,n,a,lda,ipiv,info)
    call zgetri(n,a,lda,ipiv,work,lwork,info)
    write(6,*)'info=',info
90  format(4d14.7)
    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_complex





  subroutine diag(a,b,n,lda)
    use nrtype, only : dbl,i4b
    !    use control
    use open_information 
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(n,n) :: B
    real(kind=dbl),dimension(:,:),allocatable ::vl,vr
    real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi,w
    counter2=counter2+1
    uplo='u'
    jobz='v'
    jobvl='n'
    jobvr='v'
    lda=n
    ldb=lda
    ldvr=n
    ldvl=n
    lwork=3*n-1
    itype=1
    m=lda
    ldb=lda 
    allocate(w(n)) 
    allocate(work(lwork))
    !    allocate(vr(n,n))
    !    allocate(vl(n,n))
    !    allocate(alphai(n))
    !    allocate(alphar(n))
    !    allocate(den(n))
    !    allocate(wr(n))
    !    allocate(wi(n)) 
    !    call dgegv(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,DEN,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    !    call   DGEEV( JOBVL, JOBVR, N, b, LDA, WR, WI, VL, LDVL, VR, LDVR,WORK, LWORK, INFO )
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)

    allocate(eigen_vector(n,n))
    allocate(eigen_value(n))

    eigen_vector=a
    eigen_value=w
90  format(4d14.7)
    !    write(34,*)a 
    write(23,*)w
    !    write(233,*)a
    write(6,*)'info=',info
    call cross_section(a)
    deallocate(work)
    deallocate(w)
    !    deallocate(vr)
    !    deallocate(vl)
    !    deallocate(alphai)
    !    deallocate(alphar)
    !    deallocate(den)
    deallocate(eigen_vector)
    deallocate(eigen_value)
  end subroutine diag


  subroutine diag_complex(a,b,n,lda,eigen)
    use nrtype, only : dbl,i4b 
    use control
    use open_information
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    complex*16,dimension(n,n) :: A
    complex*16,dimension(n,n) :: B
    complex*16,dimension(:,:),allocatable ::vl,vr
    complex*16,dimension(:),allocatable :: work,wr,wi,w
    real(kind=dbl),dimension(:),allocatable :: alphar,alphai,den,rwork
    complex*16,dimension(n) :: eigen 
    counter2=counter2+1
    !    counter2=2 
    uplo='u'
    jobz='v'
    jobvl='n'
    jobvr='v'
    lda=n
    ldb=lda
    ldvr=n
    ldvl=n
    lwork=2*n
    itype=1
    m=lda
    ldb=lda
    allocate(w(n))
    allocate(work(lwork))
    allocate(vr(n,n))
    allocate(vl(n,n))
    allocate(rwork(2*n))
    !allocate(alphai(n))
    !allocate(alphar(n))
    !allocate(den(n))
    !allocate(wr(n))
    !allocate(wi(n))
    !allocate(rwork(lwork))

    call ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO )

90  format(4d14.7)
    !    write(34,*)a
    !    write(35,*)w
    write(6,*)'info=',info
    eigen=w 
    ! Get the (right) eigenvectors
    A(:,:)=vr(:,:)
    deallocate(rwork)
    deallocate(work)
    deallocate(vr)
    deallocate(vl)
    deallocate(w)
    !deallocate(alphai)
    !deallocate(alphar)
    !deallocate(den)
  end subroutine diag_complex


  subroutine diag_symm(a,n,lda,eigen)
    use nrtype, only : dbl,i4b
    !    use control
    use open_information
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(:,:),allocatable ::vl,vr
    real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi,w
    real(kind=dbl),dimension(n) :: eigen
    counter2=counter2+1
    uplo='u'
    jobz='v'
    jobvl='n'
    jobvr='v'
    lda=n
    ldb=lda
    ldvr=n
    ldvl=n
    lwork=3*n-1
    itype=1
    m=lda
    ldb=lda
    allocate(w(n))
    allocate(work(lwork))
    !    allocate(vr(n,n))
    !    allocate(vl(n,n))
    !    allocate(alphai(n))
    !    allocate(alphar(n))
    !    allocate(den(n))
    !    allocate(wr(n))
    !    allocate(wi(n))
    write(6,*)'before dsygv',jobz,uplo,n,lda
    !write(370,*)a
    call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info) 
    write(6,*)'sfter dsygv info=',info
    !    allocate(eigen_vector(n,n))
    !    allocate(eigen_value(n))

    !    eigen_vector=a
    eigen=w
90  format(4d14.7)
    !    write(34,*)a
    write(236,*)w
    write(6,*)'info=',info
    !    call cross_section
    deallocate(work)
    deallocate(w)
    !    deallocate(vr)
    !    deallocate(vl)
    !    deallocate(alphai)
    !    deallocate(alphar)
    !    deallocate(den)
    !    deallocate(eigen_vector)
    !    deallocate(eigen_value)

  end subroutine diag_symm


  subroutine sort_shell(eigen_vect,open_vect,overlap)
    use nrtype
    implicit none 
    real(kind=dbl),dimension(:),INTENT(INOUT) :: eigen_vect
    real(kind=dbl),dimension(:,:),INTENT(INOUT) :: open_vect,overlap 
    real(kind=dbl),dimension(:),allocatable :: prov,provola
    integer(kind=i4b) :: i,j,inc,n,l,max
    real(kind=dbl) :: V
    n=size(eigen_vect)
    allocate(prov(n))
    allocate(provola(n))
    inc=1
    max=n
    do 
       inc=3*inc+1
       if (inc > n) exit
    end do
    do 
       inc=inc/3
       do i=inc+1,n
          v=eigen_vect(i)
          do l=1,max
             prov(l)=open_vect(l,i)
             provola(l)=overlap(l,i)
          end do
          j=i
          do 
             if (eigen_vect(j-inc) <= v) exit
             eigen_vect(j)=eigen_vect(j-inc)
             do l=1,max
                open_vect(l,j)=open_vect(l,j-inc)
                overlap(l,j)=overlap(l,j-inc)
             end do
             j=j-inc
             if(j <= inc) exit
          end do
          eigen_vect(j)=v
          do l=1,max
             open_vect(l,j)=prov(l)
             overlap(l,j)=provola(l)
          end do
       end do
       if (inc <= 1) exit
    end do
    do j=1, max
       !       write(338,*)open_vect(2,j)
       !       write(339,*)overlap(2,j)
    end do
    deallocate(prov)
    deallocate(provola)
  end subroutine sort_shell



  subroutine cross_section(H_vectors)
    use nrtype, only : dbl,i4b,nep_e,Pi,dpc
    use Brep,only : R0
    use Open_information,only : kappa,E,max_open,lmax
    use control , only : cont,molecule
    implicit none
    integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i
    real(kind=dbl) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W,Rfinal
    real(kind=dbl),dimension(:),allocatable :: J,Jp,N,Np,pre_J,pre_Jp,pre_N,pre_Np
    real(kind=dbl),dimension(:,:),allocatable :: KCC
    real(kind=dbl) :: H_vectors(max_open,max_open)
    real(kind=dbl),parameter :: ac=1.0d-14
    complex(kind=dpc) :: lambda
    !--------------------------------
    ! Subroutine that calls the ines that calculate the coulomb anr bessel
    ! functions for the matching at the sphere's boundary  
    !--------------------------------
    allocate(J(max_open))
    allocate(Jp(max_open))
    allocate(N(max_open))
    allocate(Np(max_open))
    allocate(pre_J(lmax))
    allocate(pre_Jp(lmax))
    allocate(pre_N(lmax))
    allocate(pre_Np(lmax))
    write(6,*)'molecule=',molecule
    do m=1,lmax
       xnu = m-1
       if (molecule.eq.'neutral') then
          kappa=sqrt(2.0d0*E)
          x=R0*kappa
          call sphbes(xnu, x,rj,ry,rjp,ryp)
          pre_J(m)=rj
          pre_Jp(m)=kappa*rjp
          pre_N(m)=ry
          pre_Np(m)=kappa*ryp 
       else
          x=R0*Z
          E_ryd=2.0d0*E/Z**2
          kappa=sqrt(E_ryd)
          write(6,*)'Z=',Z,R0,E_ryd,kappa


!--------------------------------------------------------------------------------------------
!!$ Coulomb functions calculation
! 1) Option1: Seaton's code
          call coulomb(cont,xnu,E_ryd,x,ac,rj,ry,rjp,ryp)
! 2) Option2: Barnett's code
!          lambda=cmplx(dble(xnu),0.d0)
!          rj=0.d0;ry=0.d0;rjp=0.d0;ryp=0.d0
!          call call_coulomb_complex(E,lambda,Z,R0,rj,ry,rjp,ryp)
!---------------------------------------------------------------------------------------------
          !*********
          Aa=1.d0
          do i=0,m-1
             Aa=Aa*(1.0d0+E_ryd*1.d0*(i**2))
          end do


          Bi=Aa/(1- exp((-2*Pi/kappa)))
!          write(6,*)Bi,E
          rj=1.d0/ sqrt(2.d0)*sqrt(Bi)*rj
          rjp=1.d0/ sqrt(2.d0)*sqrt(Bi)*rjp


          ry=1.d0/ sqrt(2.d0*Bi)*ry
          ryp=1.d0/ sqrt(2.d0*Bi)*ryp

!          write(6,*)m-1,'fl=',rj,'gl=',ry,'fjp=',rjp,'fyp=',ryp

!*COMM*** 
! Transformations to obtain Coulomb functions to match the actual wavefunction,
! not u=psi*r 

          pre_J(m) = rj/R0
          pre_Jp(m) = Z*rjp/R0-rj/(R0*R0)
          pre_N(m) = -(ry/R0)
          pre_Np(m) = -(Z*ryp/R0-ry/(R0*R0))

! *COMM*** 
! Sign is changed to the irregular coulomb function because of a different
!   convention respect to the one used by Seaton (MJ Seaton Rep. prog. phys.
!   46,167 (1983)) 


! Calculate wronskians

W=rj*ryp-ry*rjp
write(304,*)W
!*********
       endif
    end do
  channel=0
  !    princ_quant_num=lmax
    !    channel=0
    !    do k=1,princ_quant_num
    do ll=0,lmax-1
       multiplicity=2*ll+1
       do m=1,multiplicity
          channel=channel+1
          J(channel)=pre_J(ll+1)
          Jp(channel)=pre_Jp(ll+1)
          N(channel)=pre_N(ll+1)
          Np(channel)=pre_Np(ll+1)
                    write(201,*)J(channel),Jp(channel),N(channel),Np(channel),m,multiplicity,ll
       end do
    end do

    write(6,*)'after pre_J'
allocate(KCC(max_open,12))
KCC(:,:)=0.d0
    call kmatrix_perturbed(J,Jp,N,Np,Z,E,H_vectors,KCC,Rfinal)
deallocate(KCC)

    deallocate(pre_J)
    deallocate(pre_Jp)
    deallocate(pre_N)
    deallocate(pre_Np)
    deallocate(J)
    deallocate(Jp)
    deallocate(N)
    deallocate(Np)

    write(6,*)'subroutine cross_section',rj

    end subroutine cross_section


subroutine continuum(Z,R0,H,dime,eigen,KCC,ene)
  use nrtype, only : dbl,i4b,nep_e,Pi,dpc
  use Open_information,only : kappa,E,max_open,lmax
  use control , only : cont,molecule
  implicit none
  integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i,dime,kfn&
       & ,mode,ifail,ipr,jj,counter
  complex(kind=dpc) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W,beta,sig
  complex(kind=dpc),dimension(:),allocatable :: J,Jp,N,Np,pre_J,pre_Jp,pre_N,pre_Np
  real(kind=dbl),intent(inout) :: H(dime,dime),eigen(dime),R0,Z,KCC(max_open,12)
  real(kind=dbl),parameter :: ac=1.0d-14
  real(kind=dbl) :: alpha,ak,bk,ck,dk,mu,const1,ene,dummy,Rfinal,Rcalc
  complex(kind=dpc) :: lambda,W_dec,cgamma,eta
  complex(kind=dpc),parameter :: zI=(0.d0,1.d0)
  logical :: tests=.FALSE.
  external cgamma
  !--------------------------------
  ! Subroutine that calls the ones that calculate the coulomb anr bessel
  ! functions for the matching at the sphere's boundary  
  !--------------------------------
  write(6,*)'enter continuum'
  allocate(J(max_open))
  allocate(Jp(max_open))
  allocate(N(max_open))
  allocate(Np(max_open))
  allocate(pre_J(max_open))
  allocate(pre_Jp(max_open))
  allocate(pre_N(max_open))
  allocate(pre_Np(max_open))
  E=ene
  Rfinal=120.d0
  !write(6,*)'molecule=',molecule,R0,E,Z
  !-------------
  ! Tests
  !dummy=R0
  !do jj=1,20
  !R0=1.d0+(jj-1)*1.d0
  !------------------
  counter=0
  do counter=1,2
  if (counter==1) Rcalc=R0
  if (counter==2) Rcalc=Rfinal
  do m=1,max_open
     xnu = m-1
     xnu=int(sqrt(dble(m-1)))
     if (molecule.eq.'neutral') then
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! NEUTRAL PART
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (E.lt.0.d0) then
           kappa=cmplx(0.d0,sqrt(-2.0d0*E))
        else 
           kappa =cmplx(sqrt(2.d0*E),0.d0)
        end if
!write(6,*)'kappa=',kappa,E
!stop
        x=Rfinal*kappa
        rj=0.d0;ry=0.d0;rjp=0.d0;ryp=0.d0
        lambda=(-1.d0+sqrt(1.d0+4.d0*cmplx(eigen((m)),0.d0)))/2.d0
        ! Test
        !Lambda=-.5d0+zI*1.025d0
        !             
        ! Get asymptotic  bessel functions at R=Rfinal
        call sphbes(xnu, dble(x),ak,bk,ck,dk)
        !rj=ak; ry=bk; rjp=ck; ryp=dk
!---------------------------------------------------------------
!        KCC(m,1)=(ak)*sqrt(kappa)           ! If need j
!        KCC(m,2)=(kappa*ck)*sqrt(kappa)
!        KCC(m,3)=(bk)*sqrt(kappa)
!        KCC(m,4)=(kappa*dk)*sqrt(kappa)
         KCC(m,1)=(ak*Rfinal)*sqrt(kappa)           ! If need r*j
         KCC(m,2)=(kappa*ck*Rfinal+ak)*sqrt(kappa)
         KCC(m,3)=(bk*Rfinal)*sqrt(kappa)
         KCC(m,4)=(kappa*dk*Rfinal+bk)*sqrt(kappa)
!------------------------------------------------------------------

        W=KCC(m,1)*KCC(m,4)-KCC(m,2)*KCC(m,3)
        write(28,*)W,m,'bessj',Rfinal
      
        !write(32,*)xnu, dble(x),dble(rj),dble(ry),dble(rjp),dble(ryp)
        !lambda=cmplx(xnu,0.d0)
        !-----------------------------------
        ! set parameters for wclbes
        eta=(0.d0,0.d0)
        mode=1
        ! Get bessel functions
        x=Rcalc*kappa
        if (real(lambda*(lambda+1)).gt.-.25d0) then
           kfn=1
           call call_bessel_complex(x,lambda,Z,Rcalc,rj,ry,rjp,ryp,kfn)
        else
           mu = aimag (lambda)
           kfn=2
           call call_bessel_complex(x,lambda+.5d0,Z,Rcalc,rj,ry,rjp,ryp,kfn)
!------------------------------------------------------------------
           const1=sqrt(.5d0*Pi*Rcalc) ! If need r*j
           ak=const1/ sinh(.5d0*Pi*mu)*aimag(rj)
           bk=-const1/ cosh(.5d0*Pi*mu)*dble(rj)
           ck=ak/(2.d0*Rcalc) + const1/ sinh(.5d0*Pi*mu)*aimag(rjp*kappa)
           dk=(bk/(2.d0*Rcalc) - const1/ cosh(.5d0*Pi*mu)*dble(rjp*kappa))
!const1=sqrt(.5d0*Pi/Rcalc)      ! if need j
!           ak=const1/ sinh(.5d0*Pi*mu)*aimag(rj)
!           bk=-const1/ cosh(.5d0*Pi*mu)*dble(rj)
!           ck=sqrt(Pi/2.d0)/sinh(Pi*mu/2.d0)* (-.5d0/Rcalc**3/2*aimag(rj) + sqrt(1.d0/Rcalc)*aimag(rjp*kappa))
!           dk=-sqrt(Pi/2.d0)/cosh(Pi*mu/2.d0)* (-.5d0/Rcalc**3/2*dble(rj) + sqrt(1.d0/Rcalc)*dble(rjp*kappa))
!------------------------------------------------------------------

        end if
        !write(33,*)lambda, (x),(rj),(ry),(rjp),(ryp),ifail
        !----------------------------------
        if (real(lambda*(lambda+1)).gt.-.25d0) then
!------------------------------------------------------------------
!           pre_J(m)=(rj)*sqrt(kappa)           ! If need j
!           pre_Jp(m)=(kappa*rjp)*sqrt(kappa)
!           pre_N(m)=(ry)*sqrt(kappa)
!           pre_Np(m)=(kappa*ryp)*sqrt(kappa)

           pre_J(m)=(rj*Rcalc)*sqrt(kappa)           ! If need r*j
           pre_Jp(m)=(kappa*rjp*Rcalc+rj)*sqrt(kappa)
           pre_N(m)=(ry*Rcalc)*sqrt(kappa)
           pre_Np(m)=(kappa*ryp*Rcalc+ry)*sqrt(kappa)
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           write(29,*)W,m
!------------------------------------------------------------------

!,pre_J(m),pre_Jp(m),pre_N(m),pre_Np(m)

           ! Transform to Hankel function of the first kind (exponentially decaying)
           ! and take absolute value so functions are real
           !===================================================================
           ! As a quick fix take the real part of them, it is exponentially
           ! decaying for the cases looked at, lambda=.77, 4.05, might need a
           ! closer look at some point 
           !ak=dble(pre_J(m)+zI*pre_N(m))
           !bk=dble(pre_Jp(m)+zI*pre_Np(m))
           !====================================================================
!abs(dble(pre_Jp(m)+zI*pre_Np(m)))+abs(aimag(pre_Jp(m)+zI*pre_Np(m)))
           !ck=abs(pre_J(m)-zI*pre_N(m))
           !dk=abs(pre_Jp(m)-zI*pre_Np(m))
           !pre_J(m)=ak
           !pre_Jp(m)=bk
           !pre_N(m)=ck
           !pre_Np(m)=dk
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           write(28,*)W,m,'besscomp1',Rcalc

           !-------------------------------------------------------------------------
        else
           ! Positive energy solutions are identical to CWClark 1979 PRA
           ! Negative energy need to be renormalized (there is a factor of
           ! 3.99*10^-2)
           ! Get normalization factor from Wronskian for the attractive dipole
           ! channel, then normalize the matching functions 

           W=(ak*dk-ck*bk)

           pre_J(m)=ak
!/sqrt(W)       ! If have r*j
           pre_Jp(m)=ck
!/sqrt(W)
           pre_N(m)=bk
!/sqrt(W)
           pre_Np(m)=dk
!/sqrt(W)
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           write(29,*)W,m,pre_J(m),pre_Np(m),pre_N(m),pre_Jp(m),lambda,m,E
           write(28,*)W,m,'besscomp2',Rcalc

        end if

        !---------------------------------------------------------------------
        ! Tests
        if (tests) then
           ! Wronskian
           !W=rj*ryp-ry*rjp
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           write(27,*)(ak*dk-ck*bk)
           write(30,*)W,m
           ! exponential decay
           W=Rcalc
           do i=1,20
              !              R0=.1d0+.2d0*i
              Rcalc=1.d0+(i-1)*1.d0

              x=Rcalc*kappa
              !write(32,*)x,R0,kappa,lambda
              if (real(lambda*(lambda+1)).gt.-.25d0) then
                 kfn=1
                 call call_bessel_complex(x,lambda,Z,Rcalc,rj,ry,rjp,ryp,kfn)
                 Aa=(rj*Rcalc)*sqrt(kappa)           ! If need r*j
                 Bi=(ry*Rcalc)*sqrt(kappa)
                 maxo=abs(Aa+zI*Bi)
                 !write(40+m,*)R0,dble(maxo)
              else
                 kfn=2
                 call call_bessel_complex(x,lambda+.5d0,Z,Rcalc,rj,ry,rjp,ryp,kfn)
                 !write(31,*)R0,dble(rj),aimag(rj)
              end if
           end do
           Rcalc=W
        end if ! End tests
        !--------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------

     else
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! IONIC PART (Does not work)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     endif
  end do

  channel=0
  do ll=1,max_open
     !   multiplicity=2*ll+1
     !   do m=1,multiplicity
     channel=ll
     J(channel)=pre_J(ll)
     Jp(channel)=pre_Jp(ll)
     N(channel)=pre_N(ll)
     Np(channel)=pre_Np(ll)
     !write(40,*)R0,dble(J(channel)),dble(Jp(channel))

     !   end do
  end do
  ! Stick everything in KCC
  KCC(:,5+(counter-1)*4)=dble(J(:))
  KCC(:,6+(counter-1)*4)=dble(Jp(:))
  KCC(:,7+(counter-1)*4)=dble(N(:))
  KCC(:,8+(counter-1)*4)=dble(Np(:))
  !if (counter.eq.1) write(37,*)'aa',KCC(:,1),KCC(:,2),KCC(:,3),KCC(:,4)
  !if (counter.eq.1) write(37,*)'bb',KCC(:,5),KCC(:,6),KCC(:,7),KCC(:,8)
  !if (counter.eq.2) write(37,*)'cc',KCC(:,9),KCC(:,10),KCC(:,11),KCC(:,12)
  end do
  !-------------------------------------------
  ! Check the exponentially decaying functions
  !do ll=0,lmax-1
  !   lambda=(-1.d0+sqrt(1.d0+4.d0*cmplx(eigen((ll-1)**2+1),0.d0)))/2.d0
  !   alpha=aimag(lambda)
  !   beta=Pi/2+atan(tanh(Pi*alpha)*tan(alpha*log(2/sqrt(2.d0*E))-atan(aimag(cgamma(1-zI*alpha))/ real(cgamma(1-zI*alpha))))) 
  !beta=alpha*log(2.d0/sqrt(2.d0*E))-atan(aimag(cgamma(1-zI*alpha))/ real(cgamma(1-zI*alpha)))
  !     W_dec=pre_N(ll+1)*sin(Pi/sqrt(2.d0*E))+pre_J(ll+1)*cos(Pi/sqrt(2.d0*E))
  !   W_dec=-pre_N(ll+1)*sin(beta)-pre_J(ll+1)*cos(beta)
  !(pre_N(ll+1)+zI*pre_J(ll+1))*sqrt(2.d0)/2.d0
  !   write(35,*)W_dec,ll
  !   write(34,*)aimag(pre_J(ll+1))/sqrt(aimag(pre_J(ll+1))**2+aimag(pre_N(ll+1))**2),aimag(pre_J(ll+1))**2,aimag(pre_N(ll+1))**2
  !   write(34,*)sin(Pi/sqrt(2.d0*E))
  !end do
  !-------------------------------------------
  !write(6,*)'after pre_J'

  !call kmatrix(J,Jp,N,Np,E,H,dime,eigen,R)
    call kmatrix_perturbed(KCC(:,1),KCC(:,2),KCC(:,3),KCC(:,4),Z,E,H,KCC,Rfinal)

  !------------------------------
  ! End tests
  !end do
  !R0=dummy
  !stop
  !------------------------------
  deallocate(pre_J)
  deallocate(pre_Jp)
  deallocate(pre_N)
  deallocate(pre_Np)
  deallocate(J)
  deallocate(Jp)
  deallocate(N)
  deallocate(Np)

  write(6,*)'end continuum'
end subroutine continuum


  end module Solve
