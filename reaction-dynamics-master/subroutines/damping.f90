!!!*************************************************************
! 文件/File: damping.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: damping.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

! xd is some critical point 

subroutine damping(Damp,gridpts,isize)
    implicit none
    integer*8 i,isize
    real*8 xd,dx
    real*8,dimension(isize)::Damp
    real*8,dimension(1:isize),intent(in)::gridpts

    dx = 0.01
    xd = gridpts(isize) - gridpts(1) - 3.0

    do i = 1,isize
        if ( gridpts(i) .le. xd ) then
            Damp(i) = 1.0
        else
            Damp(i) = exp(-1.0*dx*(gridpts(i)-xd)**2)
        end if 
    end do

end subroutine  damping
