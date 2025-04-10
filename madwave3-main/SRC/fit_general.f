!!!*************************************************************
! 文件/File: fit_general.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: fit_general.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

      subroutine setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                       ,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension iomdiat(nelecmax),iomatom(nelecmax)
      dimension sigdiat(nelecmax),sigatom(nelecmax)

      call setxbcpot
      return
      end

      subroutine potelebond(rhAhB,rhBAu,costet,potmat,nelec,nelecmax)
      implicit real*8(a-h,o-z)
      dimension potmat(nelecmax,nelecmax)
      dimension der(3)

      return
      end
      subroutine setxbcpot
      implicit real*8(a-h,o-z)

      write(6,"(/,40('-'),//
     & ,10x,'general potential routines '
     & ,//,15x,'  , to compilae madwave3 for all systems only once'
     & ,//,40('-'),//)")

      return
      end

 
