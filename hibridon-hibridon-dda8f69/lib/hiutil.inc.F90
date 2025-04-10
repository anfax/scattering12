!!!*************************************************************
! 文件/File: hiutil.inc.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: hiutil.inc.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

! https://stackoverflow.com/questions/37500015/a-portable-way-to-suppress-an-unused-dummy-argument-warning-in-fortran
#define UNUSED(x) if (.false.) print*,shape(x)

#define IS_EXACTLY_ZERO(x) (transfer(x,1_8) == transfer(0.d0, 1_8))