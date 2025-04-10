!!!*************************************************************
! 文件/File: vibtau.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vibtau.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VIBTAU (IS, S, D, XK, PER)
C
C     VIBTAU - compute vibrational period for u motion
C
C  Called by:
C     PLAG   - compute primitive LAG probabilities from thetas
C
C  Calls:
C     WKBINT - interpolate WKB energy levels from grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      IF (LGS(7)) THEN
         CALL WKBINT (IS, S, NSTATE, E, UL, UG, PER)
      ELSE
         PER = PI*XK*0.5D0/D
         IF (LGS(3)) PER = PER/(1.D0-TNP1/XK)
      END IF
      RETURN
      END
