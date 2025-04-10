!!!*************************************************************
! 文件/File: angle.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: angle.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

module angulardvr
implicit none

contains

    subroutine angle(theta,w,mm)
    use constants

    implicit none

    integer*8 l,m,lwork,ll,mm,inf,j,k,i
    real*8,dimension(:,:),allocatable::C
    real*8,dimension(:),allocatable::eval,work
    real*8,dimension(1:ntheta)::theta,w
    real*8 sine,dif,dif1,fm,s,t1

    lwork=64*ntheta

    allocate(eval(1:ntheta))
    allocate(C(1:ntheta,1:ntheta))
    allocate(work(1:lwork))

    do m=1,ntheta
        do l=1,ntheta
            ll=l-1+mm

            if ((l-m).eq.1)then
             dif=ll**2-mm**2
             dif1=4.d0*ll**2-1
             c(l,m)=sqrt(dif/dif1)
            else if (m-l.eq.1) then
             ll=ll+1
             dif=ll**2
             dif1=4.d0*ll**2-1
             c(l,m)=sqrt(dif/dif1)
            else
             c(l,m)=0.d0
            endif

        end do
    end do

    call dsyev('V','U',ntheta,C,ntheta,eval,work,lwork,inf)
    if(inf.ne.0)print*,'dsyev error coden angle',inf

!    print*,'resulting eigenvalues- ',eval

!    print*,'eigenvector matrix is'

    do j=1,ntheta
        theta(j)=acos(eval(j))
    end do

    fm=2.d0
    if (mm.eq.0) then
        fm=fm
    else
        do i=1,mm
            fm=(2.d0*i/(2.d0*i+1))*fm
        enddo
    endif

    do i=1,ntheta
        s=dsin(theta(i))
        s=(1.d0/s)**mm
        t1=s*dabs(c(1,i))
        w(i)=fm*t1*t1
    enddo

    deallocate(C,eval,work)
end subroutine
end module
