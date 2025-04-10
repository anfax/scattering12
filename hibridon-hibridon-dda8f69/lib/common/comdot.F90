!!!*************************************************************
! 文件/File: comdot.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: comdot.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck comdot
#if defined(HIB_FPS) || defined(HIB_VAX) || defined(HIB_UNIX) || defined(HIB_CRAY) || defined(HIB_MAC)
data dot /'.'/
#endif
#if defined(HIB_UNIVAC) || defined(HIB_CRAY_COS)
data dot /'$'/
#endif
