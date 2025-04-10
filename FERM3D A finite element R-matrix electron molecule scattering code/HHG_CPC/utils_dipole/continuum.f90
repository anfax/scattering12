!!!*************************************************************
! 文件/File: continuum.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: continuum.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine read_input(nmax)
  use nrtype, only : dbl,i4b,nep_e,Pi,dpc
  use Open_information,only : kappa,E,max_open,lmax
  use control , only : cont,molecule,option_Coulomb
  use pass_quantum_numbers, only : Z,R0,DeltaE,NumSteps,Emin
  implicit none
  integer(kind=i4b) :: nmax
  open(unit=10,file='input_continuum',status='old',action='read')
  read(10,*)Z
  read(10,*)molecule
  read(10,*)R0
  read(10,*)Emin
  read(10,*)DeltaE
  read(10,*)NumSteps
  read(10,*)lmax
  close(10)
  nmax=lmax
end subroutine read_input

subroutine continuum_main(H,dime,ldb,eigen)
  use nrtype, only : dbl,i4b,nep_e,Pi,dpc
  use Open_information,only : kappa,E,max_open,lmax
  use control , only : cont,molecule,option_Coulomb
  use pass_quantum_numbers, only : Z,R0,DeltaE,NumSteps,Emin
  use root_finder
  implicit none
  integer(kind=i4b),intent(in) :: dime,ldb
  real(kind=dbl),intent(in) :: H(dime,ldb),eigen(ldb)
  integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i
  integer(kind=i4b) :: nE
  real(kind=dbl) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W,Ra,ene
  !,Z,R0
  !    real(kind=dbl) :: Emin,DeltaE
  real(kind=dbl),dimension(:),allocatable :: J,Jp,N,Np,pre_J,pre_Jp,pre_N,pre_Np
  real(kind=dbl),dimension(:,:),allocatable :: R
  real(kind=dbl),parameter :: ac=1.0d-14
  complex(kind=dpc) :: lambda
  ! Root finding vartiables
  integer(kind=i4b) :: n_brack,iR,nb
  real(kind=dbl),allocatable :: ee(:),root_final(:)
  real(kind=dbl) :: x1,x2
  !  complex(kind=dpc) ::  func_root
  real(kind=dbl),dimension(:,:),allocatable ::KCC
  !------------------------

  max_open=(lmax+1)**2
  allocate(R(max_open,max_open))
  allocate(KCC(max_open,max_open+1))
  !--------------------------------------------------------------
  ! Energy looping
  !do nE=1,NumSteps
  !   E=(Emin+(nE-1)*DeltaE)/27.2d0
  !   ene =E
  !   call R_matrix_construction(R,max_open,max_open,H,dime,eigen,E)
  !   call continuum(Z,R0,H,dime,eigen,R,KCC,ene)
  !end do
  !---------------------------------------------------------------
  !Root finding
  ! nc=max_open
  ! KCC=external, kmatrix (How?) complex
  ! ee(1:nc)=0.d0
  ! iR=1
  ! x1=Emin
  ! x2=Emax (0.d0)
  ! func_root=external (more or less already ok)
  ! n=num_bracketing intervals , try 2000, rest of the stuff allocate arrays
  !----------------------------------------------------------------
  n_brack=25000
  allocate(ee(max_open),root_final(n_brack))
  ee(:)=0.d0; iR=1; Ra=0.d0; nb=25000
  x1 = Emin
  x2 = DeltaE
  call wrap_root_finder(func_root,n_brack,KCC,ee,max_open,iR,Ra,x1,x2,nb,root_final,R,max_open,H,dime,eigen,Z,R0)
  deallocate(ee,root_final)
  !----------------------------------------------------------------
  deallocate(R,KCC)

end subroutine continuum_main

subroutine continuum(Z,R0,H,dime,eigen,R,KCC,ene)
  use nrtype, only : dbl,i4b,nep_e,Pi,dpc
  use Open_information,only : kappa,E,max_open,lmax
  use control , only : cont,molecule,option_Coulomb
  implicit none
  integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i,dime,kfn&
       & ,mode,ifail,ipr,jj
  complex(kind=dpc) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W,beta,sig
  complex(kind=dpc),dimension(:),allocatable :: J,Jp,N,Np,pre_J,pre_Jp,pre_N,pre_Np
  real(kind=dbl),intent(inout) :: H(dime,dime),eigen(dime),R0,Z,R(max_open,max_open),KCC(max_open,max_open+1)
  real(kind=dbl),parameter :: ac=1.0d-14
  real(kind=dbl) :: alpha,ak,bk,ck,dk,mu,const1,ene,dummy
  complex(kind=dpc) :: lambda,W_dec,cgamma,eta
  complex(kind=dpc),parameter :: zI=(0.d0,1.d0)
  logical :: tests=.FALSE.
  external cgamma

  !--------------------------------
  ! Subroutine that calls the ones that calculate the coulomb anr bessel
  ! functions for the matching at the sphere's boundary  
  !--------------------------------
  allocate(J(max_open))
  allocate(Jp(max_open))
  allocate(N(max_open))
  allocate(Np(max_open))
  allocate(pre_J(max_open))
  allocate(pre_Jp(max_open))
  allocate(pre_N(max_open))
  allocate(pre_Np(max_open))
  E=ene
  !write(6,*)'molecule=',molecule,R0,E,Z
  !-------------
  ! Tests
  !dummy=R0
  !do jj=1,20
  !R0=1.d0+(jj-1)*1.d0
  !------------------
  do m=1,max_open
     xnu = m-1
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
        x=R0*kappa
        rj=0.d0;ry=0.d0;rjp=0.d0;ryp=0.d0
        lambda=(-1.d0+sqrt(1.d0+4.d0*cmplx(eigen((m)),0.d0)))/2.d0
        ! Test
        !Lambda=-.5d0+zI*1.025d0
        !             
        ! test bessel functions
        !call sphbes(xnu, dble(x),ak,bk,ck,dk)
        !rj=ak; ry=bk; rjp=ck; ryp=dk
        !write(32,*)xnu, dble(x),dble(rj),dble(ry),dble(rjp),dble(ryp)
        !lambda=cmplx(xnu,0.d0)
        !-----------------------------------
        ! set parameters for wclbes
        eta=(0.d0,0.d0)
        mode=1
        ! Get bessel functions

        if (real(lambda*(lambda+1)).gt.-.25d0) then
           kfn=1
           call call_bessel_complex(x,lambda,Z,R0,rj,ry,rjp,ryp,kfn)
        else
           mu = aimag (lambda)
           kfn=2
           call call_bessel_complex(x,lambda+.5d0,Z,R0,rj,ry,rjp,ryp,kfn)

!           const1=sqrt(.5d0*Pi*R0) ! If need r*j

!           ak=const1/ sinh(.5d0*Pi*mu)*aimag(rj)
!           bk=-const1/ cosh(.5d0*Pi*mu)*dble(rj)
!           ck=ak/(2.d0*R0) + const1/ sinh(.5d0*Pi*mu)*aimag(rjp*kappa)
!           dk=(bk/(2.d0*R0) - const1/ cosh(.5d0*Pi*mu)*dble(rjp*kappa))

! Match to j
const1=sqrt(.5d0*Pi/R0)
           ak=const1/ sinh(.5d0*Pi*mu)*aimag(rj)
           bk=-const1/ cosh(.5d0*Pi*mu)*dble(rj)
           ck=sqrt(Pi/2.d0)/sinh(Pi*mu/2.d0)* (-.5d0/R0**3/2*aimag(rj) + sqrt(1.d0/R0)*aimag(rjp*kappa))
           dk=-sqrt(Pi/2.d0)/cosh(Pi*mu/2.d0)* (-.5d0/R0**3/2*dble(rj) + sqrt(1.d0/R0)*dble(rjp*kappa))

        end if
        !write(33,*)lambda, (x),(rj),(ry),(rjp),(ryp),ifail
        !----------------------------------
        if (real(lambda*(lambda+1)).gt.-.25d0) then
!           pre_J(m)=(rj*R0)*sqrt(kappa)           ! If need r*j
!           pre_Jp(m)=(kappa*rjp*R0+rj)*sqrt(kappa)
!           pre_N(m)=(ry*R0)*sqrt(kappa)
!           pre_Np(m)=(kappa*ryp*R0+ry)*sqrt(kappa)
           pre_J(m)=(rj)*sqrt(kappa)           ! If need j
           pre_Jp(m)=(kappa*rjp)*sqrt(kappa)
           pre_N(m)=(ry)*sqrt(kappa)
           pre_Np(m)=(kappa*ryp)*sqrt(kappa)

           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           !write(29,*)W,m

           ! Transform to Hankel function of the first kind (exponentially decaying)
           ! and take absolute value so functions are real
           !===================================================================
           ! As a quick fix take the real part of them, it is exponentially
           ! decaying for the cases looked at, lambda=.77, 4.05, might need a
           ! closer look at some point 
           ak=dble(pre_J(m)+zI*pre_N(m))
           bk=dble(pre_Jp(m)+zI*pre_Np(m))
           !====================================================================
!abs(dble(pre_Jp(m)+zI*pre_Np(m)))+abs(aimag(pre_Jp(m)+zI*pre_Np(m)))
           ck=abs(pre_J(m)-zI*pre_N(m))
           dk=abs(pre_Jp(m)-zI*pre_Np(m))
           pre_J(m)=ak
           pre_Jp(m)=bk
           pre_N(m)=ck
           pre_Np(m)=dk
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           !write(28,*)W,m

           !-------------------------------------------------------------------------
        else
           ! Positive energy solutions are identical to CWClark 1979 PRA
           ! Negative energy need to be renormalized (there is a factor of
           ! 3.99*10^-2)
           ! Get normalization factor from Wronskian for the attractive dipole
           ! channel, then normalize the matching functions 

           W=(ak*dk-ck*bk)

           pre_J(m)=ak/sqrt(W)       ! If have r*j
           pre_Jp(m)=ck/sqrt(W)
           pre_N(m)=bk/sqrt(W)
           pre_Np(m)=dk/sqrt(W)
           W=pre_J(m)*pre_Np(m)-pre_N(m)*pre_Jp(m)
           !write(29,*)W,m

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
           W=R0
           do i=1,20
              !              R0=.1d0+.2d0*i
              R0=1.d0+(i-1)*1.d0

              x=R0*kappa
              !write(32,*)x,R0,kappa,lambda
              if (real(lambda*(lambda+1)).gt.-.25d0) then
                 kfn=1
                 call call_bessel_complex(x,lambda,Z,R0,rj,ry,rjp,ryp,kfn)
                 Aa=(rj*R0)*sqrt(kappa)           ! If need r*j
                 Bi=(ry*R0)*sqrt(kappa)
                 maxo=abs(Aa+zI*Bi)
                 !write(40+m,*)R0,dble(maxo)
              else
                 kfn=2
                 call call_bessel_complex(x,lambda+.5d0,Z,R0,rj,ry,rjp,ryp,kfn)
                 !write(31,*)R0,dble(rj),aimag(rj)
              end if
           end do
           R0=W
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
  KCC(:,1)=dble(J(:))
  KCC(:,2)=dble(Jp(:))
  !write(37,*)KCC(:,1),KCC(:,2)
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

  !write(6,*)'subroutine cross_section',rj
end subroutine continuum



subroutine invert_real(a,n,lda)
  use control
  use nrtype, only : i4b,dbl,dpc
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
90 format(4d14.7)
  deallocate(ipiv)
  deallocate(work)
  !write(6,*)'after invert real'
end subroutine invert_real

subroutine write_simm(K,max_open,number)
  use nrtype,only : i4b,dbl,dpc,Pi
  implicit none
  integer(kind=i4b),intent(in) :: max_open,number
  real(kind=dbl), intent(in) :: K(max_open,max_open)
  integer(kind=i4b) :: p,q

  do p=1,max_open
     do q=1,max_open
        write(number,*)K(p,q),p,q
     end do
  end do

end subroutine write_simm


subroutine R_matrix_construction(R,max_open,lda,H,dime,eigen,Etot)
  use nrtype, only : i4b,dbl,dpc
  implicit none
  integer(kind=i4b) :: max_open,lda,dime
  integer(kind=i4b),save :: i,j,nchan,k,nocsf,np1,np2,kl,ir,ngeom
  integer(kind=i4b) :: counter = 0
  real(kind=dbl),intent(inout) :: R(max_open,max_open),H(dime,dime),eigen(dime),Etot
  real(kind=dbl),allocatable,save :: Rmat_UK(:,:),wamp(:,:,:),epole(:,:),rtemp(:)
  real(kind=dbl) :: sum,dummy,E_scale,Etot_in
  real(kind=dbl),save :: e0
  real(kind=dbl),parameter :: zero=0.d0
  character(10) :: option = "E-indep"
E_scale=-410.0428252280d0
Etot_in=Etot+E_scale
  call R_matrix_construction_Ap(R,max_open,lda,H,dime,eigen,Etot_in)
  call R_matrix_construction_Adp(R,max_open,lda,H,dime,eigen,Etot_in)

!--------------------------------
! Renormalization of the R-matrix
!R=R*100
!--------------------------------
!do i=1,max_open
!do j=1,max_open
!write(34,*)R(i,j),i,j
!end do
!end do

  ! Contraction
  !write(6,*)'do construction'
  R=matmul(R,H)
  R=matmul(transpose(H),R)
!-------------------
! Testing purposes
!R(:,:)=0.d0
!do j=1,max_open
!do i=1,max_open
!R(i,i)=1.d0
!write(36,*)R(i,j),i,j
!end do
!end do
!-------------------
  !write(36,*)R(:,:)
  !write(6,*)'end R-matrix construction'

end subroutine R_matrix_construction

subroutine R_matrix_construction_Ap(R,max_open,lda,H,dime,eigen,Etot)
  use nrtype, only : i4b,dbl,dpc
  implicit none
  integer(kind=i4b) :: max_open,lda,dime
  integer(kind=i4b),save :: i,j,nchan,k,nocsf,np1,np2,kl,ir,ngeom
  integer(kind=i4b),dimension(:),allocatable :: ichan
  integer(kind=i4b) :: counter = 0
  real(kind=dbl),intent(inout) :: R(max_open,max_open),H(dime,dime),eigen(dime),Etot
  real(kind=dbl),allocatable,save :: Rmat_UK(:,:),wamp(:,:,:),epole(:,:),rtemp(:)
  real(kind=dbl) :: sum,dummy
  real(kind=dbl),save :: e0
  real(kind=dbl),parameter :: zero=0.d0
  character(10) :: option = "E-indep"
  ! Read in UK R-matrix data
  nchan=15
  ir=1
  allocate(Rmat_UK(nchan,nchan),ichan(max_open))
  if (option.eq."E-dep") then
     ! Read E-dependent stuff
     do i=1,nchan
        do j=1,nchan
           read(101,*)Rmat_UK(i,j)
        end do
     end do
  else
     if (counter.eq.0) then
        read(108,*)dummy,e0
        read(104,*)
        read(104,*)ngeom,nchan,np1,np2,nocsf
        !write(6,*)ngeom,nchan,np1,np2,nocsf
        allocate(epole(np1:np2,ngeom),wamp(nchan,np1:np2,ngeom))
        allocate(rtemp(NCHAN*(NCHAN+1)*NGEOM/2))
        do I=1,NCHAN
           do  J=1,I
              K = K+1
              SUM = ZERO
              do  KL=NP1,NP2
                 if (i.eq.j) then
                    read(105,*)wamp(i,kl,ir)
                    if (i.eq.1) read(106,*) epole(kl,ir)
                 end if
              end do
           end do
        end do
     end if
     counter=counter+1
     do  K=1,NCHAN*(NCHAN+1)*NGEOM/2
        RTEMP(K) = ZERO
     end do
     !---------------------------------------------------
     ! if need to scale for the energy of the target (e0)
     !etot=etot+e0
     !---------------------------------------------------
     K=0
     DO I=1,NCHAN
        do J=1,I
           K = K+1
           SUM = ZERO
           do KL=NP1,NP2
              SUM = SUM+WAMP(I,KL,IR)*WAMP(J,KL,IR)/(EPOLE(KL,IR)-ETOT)
           end do
           RTEMP(K) = RTEMP(K) + SUM
        end do
     end do
     !write(98,*)etot
     !write(99,*)rtemp(:)
  end if
  ! Which channels present
  ! l m
  ! 1 0
  ! 3 0
  ! 3 2
  call squarm(nchan,1,rtemp,Rmat_UK)
  ! Stick them in full matrix
  R(:,:)=0.d0
  if (nchan.lt.max_open) then
ichan(:)=0
     ichan(1)=1
     ichan(2)=2
     ichan(4)=3
     ichan(5)=5
     ichan(7)=4
     ichan(9)=6
     ichan(10)=9
     ichan(12)=7
     ichan(14)=8
     ichan(16)=10
     ichan(17)=15
     ichan(19)=12
     ichan(21)=11
     ichan(23)=14
     ichan(25)=13
     do i=1,max_open
       do j=1,max_open
if ((ichan(i).le.0).or.(ichan(j).le.0)) cycle
         R(i,j)=Rmat_UK(ichan(i),ichan(j))
       end do
     end do     
  else
     R(1,1)=Rmat_UK(1,1)
     write(44,*)R(1,1),etot
  end if
  !write(34,*)Rmat_UK(1,1),Rmat_UK(2,2),Rmat_UK(3,3),Rmat_UK(1,2),Rmat_UK(1,3),Rmat_UK(2,3)
  ! Contraction
  !write(6,*)'do construction'
  !  R=matmul(transpose(H),R)
  !  R=matmul(R,H)
  deallocate(Rmat_UK,ichan)
  !write(36,*)R(:,:)
  !write(6,*)'end R-matrix construction'

end subroutine R_matrix_construction_Ap

subroutine R_matrix_construction_Adp(R,max_open,lda,H,dime,eigen,Etot)
  use nrtype, only : i4b,dbl,dpc
  implicit none
  integer(kind=i4b) :: max_open,lda,dime
  integer(kind=i4b),save :: i,j,nchan,k,nocsf,np1,np2,kl,ir,ngeom
  integer(kind=i4b) :: counter = 0
  integer(kind=i4b),dimension(:),allocatable :: ichan
  real(kind=dbl),intent(inout) :: R(max_open,max_open),H(dime,dime),eigen(dime),Etot
  real(kind=dbl),allocatable,save :: Rmat_UK(:,:),wamp(:,:,:),epole(:,:),rtemp(:)
  real(kind=dbl) :: sum,dummy
  real(kind=dbl),save :: e0
  real(kind=dbl),parameter :: zero=0.d0
  character(10) :: option = "E-indep"
  ! Read in UK R-matrix data
  ! write(6,*)'enter R_mat',counter,max_open
  nchan=10
  ir=1
  allocate(Rmat_UK(nchan,nchan),ichan(max_open))
  if (option.eq."E-dep") then
     ! Read E-dependent stuff
     do i=1,nchan
        do j=1,nchan
           read(201,*)Rmat_UK(i,j)
        end do
     end do
  else
     if (counter.eq.0) then
        read(208,*)dummy,e0
        read(204,*)
        read(204,*)ngeom,nchan,np1,np2,nocsf
        !write(6,*)ngeom,nchan,np1,np2,nocsf
        allocate(epole(np1:np2,ngeom),wamp(nchan,np1:np2,ngeom))
        allocate(rtemp(NCHAN*(NCHAN+1)*NGEOM/2))
        do I=1,NCHAN
           do  J=1,I
              K = K+1
              SUM = ZERO
              do  KL=NP1,NP2
                 if (i.eq.j) then
                    read(205,*)wamp(i,kl,ir)
                    if (i.eq.1) read(206,*) epole(kl,ir)
                 end if
              end do
           end do
        end do
     end if
     counter=counter+1
     do  K=1,NCHAN*(NCHAN+1)*NGEOM/2
        RTEMP(K) = ZERO
     end do
     !---------------------------------------------------
     ! if need to scale for the energy of the target (e0)
     !etot=etot+e0
     !---------------------------------------------------
     K=0
     DO I=1,NCHAN
        do J=1,I
           K = K+1
           SUM = ZERO
           do KL=NP1,NP2
              SUM = SUM+WAMP(I,KL,IR)*WAMP(J,KL,IR)/(EPOLE(KL,IR)-ETOT)
           end do
           RTEMP(K) = RTEMP(K) + SUM
        end do
     end do
     !write(98,*)etot
     !write(99,*)rtemp(:)
  end if
  ! Which channels present
  ! l m
  ! 1 0
  ! 3 0
  ! 3 2
  call squarm(nchan,1,rtemp,Rmat_UK)
  ! Stick them in full matrix
  !R(:,:)=0.d0
  if (nchan.lt.max_open) then
ichan(:)=0
     ichan(3)=1
     ichan(6)=2
     ichan(8)=3
     ichan(11)=5
     ichan(13)=4
     ichan(15)=6
     ichan(18)=9
     ichan(20)=7
     ichan(22)=10
     ichan(24)=8
     do i=1,max_open
       do j=1,max_open
if ((ichan(i).le.0).or.(ichan(j).le.0)) cycle
         R(i,j)=Rmat_UK(ichan(i),ichan(j))
       end do
     end do     

  else
     R(1,1)=Rmat_UK(1,1)
     write(44,*)R(1,1),etot
  end if
  !write(34,*)Rmat_UK(1,1),Rmat_UK(2,2),Rmat_UK(3,3),Rmat_UK(1,2),Rmat_UK(1,3),Rmat_UK(2,3)
  ! Contraction
  !write(6,*)'do construction'
  !  R=matmul(transpose(H),R)
  !  R=matmul(R,H)
  deallocate(Rmat_UK,ichan)
  !write(36,*)R(:,:)
  !write(6,*)'end R-matrix construction'

end subroutine R_matrix_construction_Adp

SUBROUTINE SQUARM(NDIM,NMAT,TRIANG,SQUARE)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !
  !***********************************************************************
  !
  !     SQUARM PUTS LOWER TRIANGLES BACK INTO SQUARE MATRICES)
  !
  !***********************************************************************
  !
  DIMENSION SQUARE(NDIM,NDIM,NMAT),TRIANG(*)
  !
  K = 0
  DO  L=1,NMAT
     DO  I=1,NDIM
        DO  J=1,I
           K = K+1
           SQUARE(I,J,L) = TRIANG(K)
           SQUARE(J,I,L) = TRIANG(K)
        end do
     end do
  end do
  !
END subroutine squarm

