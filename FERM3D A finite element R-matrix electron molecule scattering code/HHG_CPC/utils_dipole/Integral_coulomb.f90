!!!*************************************************************
! 文件/File: Integral_coulomb.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: Integral_coulomb.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************



module pass_quantum_numbers
  use nrtype
  integer(kind=dbl) :: n,l,m,np,lp,mp
  real(kind=dbl) :: Ene,charge
  complex(kind=dpc) :: lambda,lambdap
end module pass_quantum_numbers

subroutine  Integral_Coulomb(max_open,ss,sc,cs,cc,Z,energy,R0,J_c,Jp_c,N_c,Np_c)
  use pass_quantum_numbers,only : Ene,charge,lambda,lambdap,l,lp,m,mp
  use nrtype,only : dbl,dpc,i4b,Pi
  implicit none
  integer(kind=i4b),intent(in) :: max_open
  real(kind=dbl),intent(in) :: Z,energy,R0
  real(kind=dbl), intent(out) :: ss(max_open,max_open),sc(max_open,max_open),cs(max_open,max_open),cc(max_open,max_open),&
               & J_c(max_open),Jp_c(max_open),N_c(max_open),Np_c(max_open)
  integer(kind =i4b) :: n,nmax,dime,i,j,nn,ngp,np,lmax,mm,mmp,chan,chanp
  real(kind=dbl),allocatable :: eigen(:),H(:,:),weigh(:),xabsc(:)
  real(kind=dbl),parameter :: Debye_to_au=2.541765d0
  real(kind=dbl) :: dipole,mu,r_lo,r_hi,s_r,func_r,theta_integr,norm,s_ss,s_sc,s_cc,&
       & fc_ret,gc_ret,fcp_ret,gcp_ret,fc,gc,fcp,gcp,func_theta,theta_lo,theta_hi, &
       & func_cs,s_cs,E_ryd,Aa,Bi,kappa
  external :: func_ss,func_sc,func_cc,norm,func_theta,func_cs
! Inputs - initializations
  Ene=energy
  charge=Z
  lambda=(0.d0,0.d0)
  dipole=-4.26d0/Debye_to_au/2.d0
  r_lo=R0
  r_hi=10000.d0
  theta_lo=0.d0
  theta_hi=Pi
  ngp=1000
  E_ryd=2.d0*Ene/charge**2
  kappa=Sqrt(2.d0*Ene)

ss(:,:)=0.d0
sc(:,:)=0.d0
cs(:,:)=0.d0
cc(:,:)=0.d0
J_c(:)=0.d0
N_c(:)=0.d0
Jp_c(:)=0.d0
Np_c(:)=0.d0


!--------------------------
  write(6,*)'dipole moment=',dipole
  dime=max_open
  lmax=int(sqrt(dble(max_open)))-1
  write(6,*)'dimension=',dime,lmax
  allocate(eigen(dime))
  allocate(H(dime,dime))
  allocate(weigh(ngp))
  allocate(xabsc(ngp))
  call gauleg(ngp,xabsc,weigh)
  chan=0
  do l=0,lmax
     do mm=0,2*l
        chan=chan+1
        if (mod(mm,2).ne.0) then;m=-(mm+1)/2;else;m=(mm)/2;end if
           chanp=0
           do lp=0,lmax
              do mmp=0,2*lp
                 chanp=chanp+1
write(30,*)lp,mmp,chanp
                 if (mod(mmp,2).ne.0) then;mp=-(mmp+1)/2;else;mp=(mmp)/2;end if
                    lambda=cmplx(dble(l),0.d0)
                    lambdap=cmplx(dble(lp),0.d0)
                    if (m.ne.mp) then
                       theta_integr=0.d0
                    else if ((lp.le.l-2).or.(lp.ge.l+2)) then
                       theta_integr=0.d0
                    else if ((lp.eq.l-1).or.(lp.eq.l+1)) then
                       !call qgaus(func_theta,theta_lo,theta_hi,theta_integr,xabsc,weigh,ngp)
                       call qgaus_1d(func_theta,theta_lo,theta_hi,theta_integr)
                       !if (lp.eq.l-1) then
                       !theta_integr=2.d0*(l)/(2.d0*l-1.d0)/(2.d0*l+1.d0)*norm(l,m)*norm(lp,m)
                       !else if (lp.eq.l+1) then
                       !theta_integr=2.d0*(l+1)/(2.d0*l+1.d0)/(2.d0*l+3.d0)*norm(l,m)*norm(lp,m)
                       !end if
                       call qgaus(func_ss,r_lo,r_hi,s_ss,xabsc,weigh,ngp)
                       !write(20,*)s_ss,real(lambda),real(lambdap),m,mp,theta_integr
                       call qgaus(func_sc,r_lo,r_hi,s_sc,xabsc,weigh,ngp)
                       !write(21,*)s_sc,real(lambda),real(lambdap),m,mp,theta_integr
                       call qgaus(func_cc,r_lo,r_hi,s_cc,xabsc,weigh,ngp)
                       !write(22,*)s_cc,real(lambda),real(lambdap),m,mp,theta_integr
                       call qgaus(func_cs,r_lo,r_hi,s_cs,xabsc,weigh,ngp)
                       !write(23,*)s_cs,real(lambda),real(lambdap),m,mp,theta_integr
                       write(24,*)s_ss,real(lambda),real(lambdap),m,mp,theta_integr,l,lp
                       call call_coulomb_complex(Ene,lambdap,charge,R0,fc_ret,gc_ret,fcp_ret,gcp_ret)
                       write(20,*)s_ss,real(lambda),real(lambdap),m,mp,theta_integr,chan,chanp,mm,mmp
                       ss(chan,chanp)=s_ss*theta_integr*dipole
                       sc(chan,chanp)=s_sc*theta_integr*dipole
                       cs(chan,chanp)=s_cs*theta_integr*dipole
                       cc(chan,chanp)=s_cc*theta_integr*dipole
                    end if
                 end do
              end do
           end do
        end do

        chan=0
        chanp=0
        do l=0,lmax
          lambda=cmplx(dble(l),0.d0)
           do mm=0,2*l
              chan=chan+1
              if (mod(mm,2).ne.0) then;m=-(mm+1)/2;else;m=(mm)/2;end if
                 call call_coulomb_complex(Ene,lambda,charge,R0,fc_ret,gc_ret,fcp_ret,gcp_ret)
                 fc=fc_ret
                 gc=gc_ret
                 fcp=fcp_ret
                 gcp=gcp_ret
                 do chanp=1,max_open
                    fc=fc+  Pi*(gc_ret *ss(chan,chanp)-fc_ret *cs(chan,chanp))
                    gc=gc+  Pi*(gc_ret *sc(chan,chanp)-fc_ret *cc(chan,chanp))
                    fcp=fcp+Pi*(gcp_ret*ss(chan,chanp)-fcp_ret*cs(chan,chanp))
                    gcp=gcp+Pi*(gcp_ret*sc(chan,chanp)-fcp_ret*cc(chan,chanp))
                 end do
                 !----------------------------------------------

write(307,*)Ene,lambda,charge,R0,fc_ret,gc_ret,fcp_ret,gcp_ret,lambda
write(308,*)Ene,lambda,charge,R0,fc,gc,fcp,gcp,lambda

                 Aa=1.d0
                 do i=0,l
                    Aa=Aa*(1.0d0+E_ryd*1.d0*(i**2))
                 end do

                 Bi=Aa/(1- exp((-2*Pi/kappa)))
                 fc=1.d0/ sqrt(2.d0)*sqrt(Bi)*fc
                 fcp=1.d0/ sqrt(2.d0)*sqrt(Bi)*fcp

                 gc=1.d0/ sqrt(2.d0*Bi)*gc
                 gcp=1.d0/ sqrt(2.d0*Bi)*gcp

                 J_c(chan) = fc/R0
                 Jp_c(chan) = Z*fcp/R0-fc/(R0*R0)
                 N_c(chan) = -(gc/R0)
                 Np_c(chan) = -(Z*gcp/R0-gc/(R0*R0))
write(202,*)J_c(chan),Jp_c(chan),N_c(chan),Np_c(chan),l,chan
                 !----------------------------------------------
           end do
        end do
           deallocate(weigh)
           deallocate(xabsc)
           deallocate(H)
           deallocate(eigen)
         end subroutine Integral_Coulomb

         function func_ss(rr)
           use nrtype, only : dbl,i4b,Pi
           use pass_quantum_numbers
           implicit none
           real(kind=dbl) :: func_ss,lag_pol(0:n+l),lag_polp(0:np+lp),factrl,rescaled_variable,norm,normp,&
                & rescaled_variablep,fc_ret,gc_ret,fcp_ret,gcp_ret,xx,fc_ret1,fc_ret2
           real(kind=dbl),intent(in) :: rr
           call call_coulomb_complex(Ene,lambda,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret1=fc_ret
           call call_coulomb_complex(Ene,lambdap,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret2=fc_ret
           func_ss=fc_ret1*fc_ret2/rr**2
           !xx=Sqrt(2.d0*Ene)*rr
           !write(70,*)fc_ret,rr
           !xx/2.d0
           !write(71,*)func_r,rr
           !write(72,*)fc_ret*gc_ret,rr
           !write(73,*)gc_ret*gc_ret,rr
           !write(74,*)gc_ret,rr
           !xx/2.d0
         end function func_ss

         function func_sc(rr)
           use nrtype, only : dbl,i4b,Pi
           use pass_quantum_numbers
           implicit none
           real(kind=dbl) :: func_sc,lag_pol(0:n+l),lag_polp(0:np+lp),factrl,rescaled_variable,norm,normp,&
                & rescaled_variablep,fc_ret,gc_ret,fcp_ret,gcp_ret,xx,fc_ret1,fc_ret2
           real(kind=dbl),intent(in) :: rr
           call call_coulomb_complex(Ene,lambda,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret1=fc_ret
           call call_coulomb_complex(Ene,lambdap,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret2=gc_ret
           func_sc=fc_ret1*fc_ret2/rr**2
         end function func_sc

         function func_cs(rr)
           use nrtype, only : dbl,i4b,Pi
           use pass_quantum_numbers
           implicit none
           real(kind=dbl) :: func_cs,lag_pol(0:n+l),lag_polp(0:np+lp),factrl,rescaled_variable,norm,normp,&
                & rescaled_variablep,fc_ret,gc_ret,fcp_ret,gcp_ret,xx,fc_ret1,fc_ret2
           real(kind=dbl),intent(in) :: rr
           call call_coulomb_complex(Ene,lambda,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret1=gc_ret
           call call_coulomb_complex(Ene,lambdap,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret2=fc_ret
           func_cs=fc_ret1*fc_ret2/rr**2
         end function func_cs


         function func_cc(rr)
           use nrtype, only : dbl,i4b,Pi
           use pass_quantum_numbers
           implicit none
           real(kind=dbl) :: func_cc,lag_pol(0:n+l),lag_polp(0:np+lp),factrl,rescaled_variable,norm,normp,&
                & rescaled_variablep,fc_ret,gc_ret,fcp_ret,gcp_ret,xx,fc_ret1,fc_ret2
           real(kind=dbl),intent(in) :: rr
           call call_coulomb_complex(Ene,lambda,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret1=gc_ret
           call call_coulomb_complex(Ene,lambdap,charge,rr,fc_ret,gc_ret,fcp_ret,gcp_ret)
           fc_ret2=gc_ret
           func_cc=fc_ret1*fc_ret2/rr**2
         end function func_cc

         function func_theta(theta)
           use nrtype, only : dbl,i4b,Pi
           use pass_quantum_numbers
           implicit none
           real(kind=dbl) :: func_theta,pleg,factrl
           real(kind=dbl),intent(in) :: theta
           external :: factrl,pleg
           func_theta=pleg(l,m,theta)*pleg(lp,mp,theta)*sin(theta)*2.d0*Pi*cos(theta)
         end function func_theta

         subroutine diag_symm(a,n,lda,eigen)
           use nrtype, only : dbl,i4b
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
           write(6,*)'before dsygv',jobz,uplo,n,lda
           call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
           write(6,*)'sfter dsygv info=',info

           eigen=w
90         format(4d14.7)
           write(23,*)w
           write(6,*)'info=',info
           deallocate(work)
           deallocate(w)

         end subroutine diag_symm

         subroutine diag_complex(a,n,lda,eigen)
           use nrtype, only : dbl,i4b
           implicit none
           integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
           character jobvl,jobvr,jobz,uplo
           complex*16,dimension(n,n) :: A
           complex*16,dimension(:,:),allocatable ::vl,vr
           complex*16,dimension(:),allocatable :: work,wr,wi,w
           real(kind=dbl),dimension(:),allocatable :: alphar,alphai,den,rwork
           complex*16,dimension(n) :: eigen
           counter2=counter2+1
           uplo='u'
           jobz='v'
           jobvl='n'
           jobvr='v'
           lda=max(1,n)
           ldb=lda
           ldvr=n
           ldvl=n
           lwork=max(1,2*n)
           itype=1
           m=lda
           ldb=lda
           allocate(w(n))
           allocate(work(lwork))
           allocate(vr(n,n))
           allocate(vl(n,n))
           allocate(alphai(n))
           allocate(alphar(n))
           allocate(den(n))
           allocate(wr(n))
           allocate(wi(n))
           allocate(rwork(lwork))
           call ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO )
           if (info.ne.0) then
              write(6,*)'info=',info
           end if
           eigen=w
           deallocate(work)
           deallocate(vr)
           deallocate(vl)
           deallocate(alphai)
           deallocate(alphar)
           deallocate(den)
         end subroutine diag_complex


         SUBROUTINE matinv_real(mattmp,nen2,lda)
           use nrtype, only : dbl
           IMPLICIT NONE

           INTEGER,intent(IN) :: nen2,lda
           real(kind=dbl),intent(INOUT) :: mattmp(nen2,nen2)

           INTEGER :: i,j,info,lwork
           INTEGER,ALLOCATABLE :: ipiv(:)
           real(kind=dbl),ALLOCATABLE :: work(:),mattmp2(:,:)

           lwork=3*nen2
           allocate(ipiv(nen2),work(lwork))


           call dgetrf(nen2,nen2,mattmp,nen2,ipiv,info)
           call dgetri(nen2,mattmp,nen2,ipiv,work,lwork,info)

           deallocate(ipiv,work)
         END SUBROUTINE matinv_real


         subroutine modify_coulomb(R0,Z,E,max_open,lmax,J,Jp,N,Np)
           use nrtype, only : dbl,i4b,nep_e,Pi,dpc
           implicit none
           integer(kind=i4b),intent(in) :: lmax,max_open
           real(kind=dbl),intent(in) :: R0,Z,E
           real(kind=dbl),intent(inout) :: J(max_open),Jp(max_open),N(Max_open),Np(max_open)
           integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i
           real(kind=dbl) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W,kappa
           real(kind=dbl),dimension(:),allocatable :: pre_J,pre_Jp,pre_N,pre_Np
           real(kind=dbl),parameter :: ac=1.0d-14
           complex(kind=dpc) :: lambda
           character*15 :: molecule='neutral'
           logical :: cont=.TRUE.
           !--------------------------------
           ! Subroutine that calls the ines that calculate the coulomb anr bessel
           ! functions for the matching at the sphere's boundary
           !--------------------------------
           allocate(pre_J(lmax))
           allocate(pre_Jp(lmax))
           allocate(pre_N(lmax))
           allocate(pre_Np(lmax))
           do m=1,lmax
              xnu = m-1
                 x=R0*Z
                 E_ryd=2.0d0*E/Z**2
                 kappa=sqrt(E_ryd)
                 lambda=cmplx(dble(xnu),0.d0)
                 call call_coulomb_complex(E,lambda,Z,R0,rj,ry,rjp,ryp)
                 Aa=1.d0
                 do i=0,m-1
                    Aa=Aa*(1.0d0+E_ryd*1.d0*(i**2))
                 end do

                 Bi=Aa/(1- exp((-2*Pi/kappa)))
                 rj=1.d0/ sqrt(2.d0)*sqrt(Bi)*rj
                 rjp=1.d0/ sqrt(2.d0)*sqrt(Bi)*rjp

                 ry=1.d0/ sqrt(2.d0*Bi)*ry
                 ryp=1.d0/ sqrt(2.d0*Bi)*ryp
                 pre_J(m) = rj/R0
                 pre_Jp(m) = Z*rjp/R0-rj/(R0*R0)
                 pre_N(m) = -(ry/R0)
                 pre_Np(m) = -(Z*ryp/R0-ry/(R0*R0))

                 W=rj*ryp-ry*rjp
                 write(304,*)W
           end do
           channel=0
           do ll=0,lmax-1
              multiplicity=2*ll+1
              do m=1,multiplicity
                 channel=channel+1
                 J(channel)=pre_J(ll+1)
                 Jp(channel)=pre_Jp(ll+1)
                 N(channel)=pre_N(ll+1)
                 Np(channel)=pre_Np(ll+1)
              end do
           end do
           deallocate(pre_J)
           deallocate(pre_Jp)
           deallocate(pre_N)
           deallocate(pre_Np)

           write(6,*)'subroutine modify_coulomb',rj

           end subroutine modify_coulomb
