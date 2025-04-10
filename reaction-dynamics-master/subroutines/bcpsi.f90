!!!*************************************************************
! 文件/File: bcpsi.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: bcpsi.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

    subroutine bcpsi(n,x,delr1,ev,H)

    implicit none

    real*8 m,sum1,delr1,delx
    real*8,parameter::pi=3.1416,xe=1.40201
    integer*8 n,i,i_,l,inf
    real*8,dimension(:,:),allocatable::T,eve
    real*8,dimension(:),allocatable::V,work
    real*8,dimension(1:n,1:n)::H
    real*8,dimension(1:n)::ev,x
    real*8,dimension(1:3)::rad
    real*8 amu,bmass,cmass,rm

    !mass
    bmass=1.007
    cmass=1.007

    amu = 1822.89d0
    rm = amu*bmass*cmass/(bmass+cmass)

    l=3*n+3
    delx=(x(n)-x(1))/(n-1)

    allocate(T(1:N,1:N))
    allocate(eve(1:N,1:N))
    allocate(V(1:N))
    allocate(work(1:l))

    rad(1)=100.0

    open (unit=1,file='bcpot.out')
    open (unit=3,file='bchamil.out')

    do i=1,n
        rad(2)=x(i)
        rad(3)=rad(2)+rad(1)
        call bkmp2(rad,V(i))
        v(i) = v(i) + 0.17449577081d0
        write(1,*),x(i),v(i)
    end do

      do i=1,n
        do i_=1,n
            if (i.eq.i_) then
                T(i,i_)=(((pi**2/3.0)-(1.0/(2.0*(i**2))))/(2.0*rm*(delx**2)))
                H(i,i_)=T(i,i_)+V(i)
            else
                T(i,i_)=((-1)**(i-i_))*((2.0/((i-i_)**2))-(2.0/((i+i_)**2)))
                T(i,i_)=T(i,i_)/(2.0*rm*(delx**2))
                H(i,i_)=T(i,i_)
            end if
        end do
        write(3,*),(H(i,i_),i_=1,n)
    end do
    close(1)
    close(3)

    call dsyev('V','U',n,H,n,ev,work,l,inf)
    print*,'dsyev error coden bcpsi',inf
    do i=1,2
     print*,ev(i)*27.2116d0
    end do


end subroutine
