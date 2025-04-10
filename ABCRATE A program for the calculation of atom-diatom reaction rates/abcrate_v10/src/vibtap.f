!!!*************************************************************
! 文件/File: vibtap.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vibtap.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VIBTAP (NF, S, D, XK, PER)
C
C     VIBTAP - compute vibrational period for u motion in product
C        channel
C
C  Called by:
C     PLAG   - compute primitive LAG probabilities from thetas
C
C  Calls:
C     AITKF2 - Aitken interpolation
C     LOCS   - locate position of S in grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(6)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      COMMON /VIBTA2/ TAUS(NSDM1)
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE NINT, IS                                                     TCA1097
      DATA NINT /5/, IS /0/
C
      IF (LGS(3).AND.LGS(7)) THEN
C interpolate period from gird of WKB values. (These are set in VSPLN2)
         IF (S .LE. SS(ISN)) THEN
            PER = TAUS(ISN)
         ELSE IF (S .GE. SS(NS2)) THEN
            PER = TAUS(NS2)
         ELSE
            CALL LOCS (IS, S)
            NP = NINT+1
            NN = NP/2 - 1
            IS = MAX(1, IS-NN)
            IF (IS+NP .GT. NS2) IS = NS2 - NP
            PER = AITKF2(S, TAUS(IS), F, SS(IS), NINT)
         END IF
C        WRITE (99, 9900) IS, ISN, NS2, S, TAUS(IS), PER
C9900    FORMAT (' IS,ISN,NS2,S,TAUS,PER=', 3I5, 1P3E15.7)
      ELSE
C compute from harmonic or Morse approximation
         PER = PI*XK*0.5D0/D
         IF (LGS(3)) PER = PER/(1.D0-DBLE(2*NF+1)/XK)                   GCL1092
      END IF
      RETURN
      END
