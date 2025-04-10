!!!*************************************************************
! 文件/File: sagarg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: sagarg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SAGARG (PS, S, D, XK, CP, DUDSM, A)
C
C     SAGARG - compute effective mass terms for adiabatic tunneling
C              calcuations
C
C  Called by:
C     PSAG   - compute SAG-type probabilites
C
C  Calls:
C     SAGCMP - compute info needed for effective mass terms
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)                                                    GCL0493
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      RPS = SQRT(PS)
C  MEPSAG
      ISAG = 1
      A(ISAG) = RPS
      CALL SAGCMP(S, D, XK, CP, AM, OM, XLM, UM)
C  EMMCPSAG
      ISAG = ISAG + 1
      CPUM = CP*UM
      T2 = 1.0D0 - CPUM
      T2 = MAX(0.0D0, T2)
      D2 = DUDSM*DUDSM
      A(ISAG) = SQRT(T2*T2+D2)*RPS
C  SCSAG
      ISAG = ISAG + 1
      T = -CPUM*(1.0D0+0.5D0*CPUM) + 0.5D0*D2
      A(ISAG) = RPS
      IF (T .LT. 0.0D0) A(ISAG) = RPS*EXP(T)
C  PASAG
      ISAG = ISAG + 1
      T = 1.0D0 - XLM
      A(ISAG) = 0.0D0
      IF (T .GT. 0.0D0) THEN
         T = T*SQRT(T)
         IF (LGS(3)) T = T/(1.D0+3.D0*CP*TNP1/(AM*XK))
         A(ISAG) = SQRT(T)*RPS
      END IF
      RETURN
      END
