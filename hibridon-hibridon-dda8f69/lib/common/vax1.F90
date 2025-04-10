!!!*************************************************************
! 文件/File: vax1.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: vax1.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck vax1
#if defined(HIB_UNIX_SUN)
      integer shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_VAX) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_NEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_DEC)
      integer xor,or,and,shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
#if defined(HIB_UNIX_IRIS)
      integer xor,or,and,shiftl,shiftr,popcnt,poppar,compl
      logical mnpar
#endif
