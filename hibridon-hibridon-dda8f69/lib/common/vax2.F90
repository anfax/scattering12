!!!*************************************************************
! 文件/File: vax2.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: vax2.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck vax2
#if defined(HIB_UNIX_AIX) || defined(HIB_VAX) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_NEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_DEC)
      xor(i,j)=ieor(i,j)
      or(i,j)=ior(i,j)
      and(i,j)=iand(i,j)
      shiftl(i,j)=ishft(i,j)
      shiftr(i,j)=ishft(i,-j)
#endif
#if defined(HIB_UNIX_IRIS)
      xor(i,j)=ieor(i,j)
      or(i,j)=ior(i,j)
      and(i,j)=iand(i,j)
      shiftl(i,j)=ishft(i,j)
      shiftr(i,j)=ishft(i,-j)
#endif
#if defined(HIB_UNIX_SUN)
      shiftl(i,j)=lshift(i,j)
      shiftr(i,j)=lrshft(i,j)
      compl(i)=not(i)
#endif
#if defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_NEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_AIX) || defined(HIB_IBM) || defined(HIB_UNIX_DEC)
      compl(i)=not(i)
#endif
#if defined(HIB_UNIX_IRIS)
      compl(i)=not(i)
#endif
#if defined(HIB_UNIX_CONVEX)
      locate(li3562,l3562,i3562,ip3562)=iisveq(l3562,li3562,i3562,ip3562)
#endif
