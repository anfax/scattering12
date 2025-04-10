!!!*************************************************************
! 文件/File: ortog.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: ortog.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


subroutine ortog (z,b,n,orto)
  use open_information
  implicit none
  integer i,j,f,n,k,l
  real*8 :: countero,proof,epsilon_custom
  real*8, dimension (:,:),allocatable :: psi
  real*8, dimension (:), allocatable :: N2,coeff
  real*8, dimension (n,n) :: z,b,orto
  write(327,*)orto
  allocate(coeff(n))
  allocate(psi(n,n))
  allocate(N2(n))
  epsilon_custom=.00001
  psi=0.0d0
  coeff=0.0d0
  N2=0.0d0
  write(326,*)orto 
  do l=1,n 
     countero=0
     do k=1,n
        if ((orto(l,k).gt.epsilon_custom).and.(k.ne.l)) then
           countero=countero+1
           N2(l)=0.0d0
           do j=1, n
              do f=1,n
                 N2(l)=N2(l)+z(j,l)*z(f,l)*b(j,f)
              end do
           end do
           write(271,*)N2(l)
           if (N2(l).lt. 0.0d0) then
              N2(l)=-N2(l)
           end if
           do j=1,n 
              psi(j,l)=z(j,l)/sqrt(N2(l))
           end do
           do j=1,n
              do f=1,n
                 coeff(l)=coeff(l)-(z(j,k)*z(f,l)*b(j,f))/sqrt(N2(l))
              end do
           end do
           proof=0.0d0
           do j=1,n
              do f=1,n
                 proof=proof+z(j,k)*b(j,f)*psi(f,l)
              end do
           end do
           proof=proof+coeff(l)
           do j=1,n  
              z(j,l)=z(j,l)/sqrt(N2(l))
              z(j,k)=z(j,k)+z(j,l)*coeff(l)
           end do
           N2(k)=0.0d0 
           do j=1, n
              do f=1,n
                 N2(k)=N2(k)+z(j,k)*z(f,k)*b(j,f)
              end do
              if (N2(k).lt.0) then
                 N2(k)=-N2(k)
              end if
              psi(j,k)=psi(j,k)/sqrt(N2(k))
           end do
           do j=1,n
              z(j,k)=z(j,k)/sqrt(N2(k))
           end do

        else if ((k.eq.n).and.(countero.eq.0)) then
           do j=1, n
              do f=1,n
                 N2(l)=N2(l)+z(j,l)*z(f,l)*b(j,f)
              end do
           end do
           if (N2(l).lt.0) then
              N2(l)=-N2(l)
           end if
           do j=1,n
              z(j,l)=z(j,l)/sqrt(N2(l))
           end do
        end if
     end do
  end do
  deallocate(coeff)
  deallocate(psi)
  deallocate(N2)
  write(314,*)proof
end subroutine ortog

