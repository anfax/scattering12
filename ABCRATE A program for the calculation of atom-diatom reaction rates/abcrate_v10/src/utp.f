!!!*************************************************************
! 文件/File: utp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: utp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      FUNCTION UTP (D, XK, XMU)
C
C     UTP    - find outer turning point in u motion
C
C  Called by:
C     DERS   - derivatives of morse turning point and zero pt. energy
C     SAGCMP - compute info needed for effective mass terms
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      IF (LGS(3)) THEN
C  Morse
         AM = SQRT(8.D0*XMU*D)/XK
         T = 1.D0 - TNP1/XK
         UTP = 1000.D0                                                  GCL1092
         IF (T .GT. 0.0D0) UTP = -LOG(1.D0-SQRT(1.D0-T*T))/AM
      ELSE
C  Harmonic
         UTP = 0.5D0*SQRT(TNP1*XK/(XMU*D))
      END IF
      RETURN
      END
