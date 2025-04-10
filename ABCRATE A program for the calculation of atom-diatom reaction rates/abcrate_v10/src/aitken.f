!!!*************************************************************
! 文件/File: aitken.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: aitken.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      SUBROUTINE AITKEN (IS, S, V, D, XK, CP)
C
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR from r.p.
C        grid
C
C  Called by:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C     LALFST - set up zeroth order path for LAG
C     LGDS   - compute velocity factor and vibrational period for LAG
C     PLAG   - compute primitive LAG probabilities from thetas
C     PSAG   - compute SAG-type probabilites
C
C  Calls:
C     AITKF2 - Aitken interpolation
C     LOCS   - locate position of s in grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION F(7)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /INTER1/ NINT
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL1096
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      DO 200 I = 1, 7                                                   GCL1096
             F(I) = 0.D0                                                GCL1096
200   CONTINUE                                                          GCL1096
C
      IF (S .LE. SS(1)) THEN
C  If S is off the grid set values to last point on grid
         V = VS(1)
         D = DS(1)
         XK = XKS(1)
         CP = CAPS(1)
      ELSEIF (S .GE. SS(NSMAX)) THEN
         V = VS(NSMAX)
         D = DS(NSMAX)
         XK = XKS(NSMAX)
         CP = CAPS(NSMAX)
      ELSE
C  Interpolate
         NP = NINT + 1
         NN = NP/2 - 1
C  Locate grid point for S on the grid
         CALL LOCS (IS, S)
         IS = MAX(1, IS-NN)
         IS = MIN(IS, NSMAX - NP)
         V = AITKF2(S, VS(IS), F, SS(IS), NINT)
         D = AITKF2(S, DS(IS), F, SS(IS), NINT)
         XK = AITKF2(S, XKS(IS), F, SS(IS), NINT)
C  check that the curvature is not zero at middle grid point; if so set
C     interpolation of curvature to zero.
         IF (CAPS(IS+NN) .NE. 0.0D0) THEN
            CP = AITKF2(S, CAPS(IS), F, SS(IS), NINT)
         ELSE
            CP = 0.0D0
         END IF
      END IF
      RETURN
      END
