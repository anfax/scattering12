!!!*************************************************************
! 文件/File: estr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: estr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      FUNCTION ESTR (IS, S, D, XK)
C
C     ESTR   - compute stretch vibrational energy levels
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     GTST   - compute free energies, CVT and ICVT rates
C     RPHSUM - summarize reaction path info
C     TABL21 - print out table of GTS info
C     THRCOR - compute threshold corrections for ICVT
C
C  Calls:
C     WKBINT - interpolate WKB energies from grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      IF (LGS(7)) THEN
         CALL WKBINT(IS, S, NSTATE, E, UL, UG, PER)
      ELSE
         E = 2.D0*TNP1*D/XK
         IF (LGS(3)) THEN
            IF (TNP1 .GE. XK) THEN
               E = D
            ELSE
               E = E*(1.D0-0.5D0*TNP1/XK)
            END IF
         END IF
      END IF
      ESTR = E
      RETURN
      END
