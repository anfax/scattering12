!!!*************************************************************
! 文件/File: parlbuf.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: parlbuf.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck parlbuf
#if defined(HIB_VAX) || defined(HIB_MAC) || defined(HIB_FPS)
parameter (lbuf=512)
#endif
#if defined(HIB_UNIX_DEC)
parameter (lbuf=1024)
#endif
