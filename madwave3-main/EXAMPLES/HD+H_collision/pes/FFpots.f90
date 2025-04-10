!!!*************************************************************
! 文件/File: FFpots.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: FFpots.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

subroutine setpotxyz
implicit none

write(*,*) 'H3pot'
write(*,*) 'No analytical derivatives yet'
endsubroutine

subroutine VH3(RAB,RBC,RAC,VTOT,DVTOT,ID)
use H3_pot
implicit none
real(8), intent(in) :: RAB,RBC,RAC
real(8), intent(out) :: VTOT
real(8), dimension(3), intent(out) :: DVTOT 
integer, intent(out) :: ID

! For derivatives:
integer :: i
real(8), dimension(3) :: R, R_tmp
real(8) :: VTOT1, VTOT2
real(8), parameter :: delta=1e-6

call NN_VH3(RAB,RBC,RAC,VTOT, DVTOT)

end subroutine
