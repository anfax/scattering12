!!!*************************************************************
! 文件/File: fumcp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: fumcp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      FUNCTION FUMCP (S,D,XK,CP,LNONAD,NFINAL)
C
C     FUMCP  - compute outer turning point for u motion (Marcus-Coltrin
C              path)
C
C  Called by:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C     LALFST - set up zeroth order path for LAG
C
C  Calls:
C     WKBINT - interpolate WKB energy levels from grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LNONAD
      DIMENSION F(6)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      COMMON /TP2CM/ ULS(NSDM1),UGS(NSDM1)
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE NINT, IS                                                     TCA1097
      DATA NINT /5/, IS /0/
C
      IF (LNONAD) THEN
         TNP1P = 2.0D0*DBLE(NFINAL) + 1.0D0                             GCL1092
      ELSE
         TNP1P = TNP1
      END IF
      SGN = SIGN(1.D0, CP)
      IF (LGS(3).AND.LGS(7)) THEN
C  WKB turning points
         IF (LNONAD) THEN
C  excited state 
            IF (S .LE. SS(1)) THEN
               UL = ULS(1)
               UG = UGS(1)
            ELSE IF (S .GE. SS(NS2)) THEN
               UL = ULS(NS2)
               UG = UGS(NS2)
            ELSE
               CALL LOCS (IS, S)
               NP = NINT+1
               NN = NP/2 - 1
               IS = MAX(1, IS-NN)
               IF (IS+NP .GT. NS2) IS = NS2 - NP
               UL = AITKF2(S, ULS(IS), F, SS(IS), NINT)
               UG = AITKF2(S, UGS(IS), F, SS(IS), NINT)
            END IF
         ELSE
C  use initial state
            CALL WKBINT(0, S, NSTATE, E, UL, UG, PER)
         END IF
      ELSE IF (LGS(3)) THEN
C  Morse approximation
         T1 = 1.D0-TNP1P/XK
         T2 = SQRT(8.D0*XMU*D)/XK
         T1 = SQRT(1.D0-T1*T1)
         UL = -LOG(1.D0+T1)/T2
         UG = -LOG(1.D0-T1)/T2
      ELSE
C  Harmonic
         UL = -0.5D0*SQRT(TNP1P*XK/(XMU*D))
         UG = -UL
      END IF
      IF (SGN .GT. 0.D0) THEN
         FUMCP = UG
      ELSE
         FUMCP = UL
      END IF
      RETURN
      END
