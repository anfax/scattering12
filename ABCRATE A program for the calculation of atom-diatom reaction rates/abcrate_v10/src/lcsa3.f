!!!*************************************************************
! 文件/File: lcsa3.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lcsa3.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LCSA3 (CSA, ALF, DS, DT2, TDV)
C
C     LCSA3  - compute cosine factors for LAG calculation
C
C  Called by:
C     LCSA2  - compute cosine factors for LAG calculation
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ALFM = 1.D0 - ALF
      DS2 = DS*DS
      T3 = SQRT(ALFM*ALFM*DS2 + ALF*(ALF*DT2 + 2.D0*ALFM*DS*TDV))
      CSA = (ALFM*DS + ALF*TDV)/T3
      RETURN
      END
