!!!*************************************************************
! 文件/File: sagcmp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: sagcmp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SAGCMP (S, D, XK, CP, AM, OM, XLM, UM)
C
C     SAGCMP - compute info needed for effective mass terms
C
C  Called by:
C     SAGARG - compute effective mass terms for adiabatic tunneling
C              calcuations
C
C  Calls:
C     WKBINT - interpolate WKB energy levels from grid
C     UTP    - find outer turning point in u motion
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      AM = SQRT(8.0D0*XMU*D)/XK
      OM = 4.0D0*D/XK
      XLM = CP*CP*TNP1/(XMU*OM)
      IF (LGS(7)) THEN
         CALL WKBINT (0, S, NSTATE, E, UL, UM, PER)
      ELSE
         UM = UTP (D, XK, XMU)
      END IF
      RETURN
      END
