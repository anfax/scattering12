!!!*************************************************************
! 文件/File: momentGrid.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: momentGrid.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

subroutine makeMomentumGrid(Kgrid,npts,length)

    implicit none
    integer*8 npts,i,m
    real*8,dimension(1:npts)::Kgrid
    real*8 length,N,pi

    pi = 3.141592653589
    N = real(npts)
    do i=1,npts
        m = i-1
        if (m.lt.N/2.0) then
            kgrid(i) = (2*pi*m)/length
        else if (m.ge.N/2.0) then
            kgrid(i) = (2*pi*(m-N))/length
        endif
    end do
    
end subroutine makeMomentumGrid