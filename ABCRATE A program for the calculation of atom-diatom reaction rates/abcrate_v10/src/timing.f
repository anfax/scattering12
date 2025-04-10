!!!*************************************************************
! 文件/File: timing.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: timing.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TIMING (TIMEF, TIMEI)
C
C     TIMING - prints out timing information, current and elasped time.
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C     PORCPU - returns cpu time
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL PORCPU(TIMEF)                                                GCL0992
      TIME = TIMEF-TIMEI
      WRITE (6, 601) TIMEF, TIME
      TIMEI = TIMEF
      RETURN
601   FORMAT(/,1X,78('*'),/,1X,T5,'Time is',F12.2,' seconds, elapsed ',
     *             'time is',F12.2,' seconds',/,1X,78('*'),/)
      END
