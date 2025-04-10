!!!*************************************************************
! 文件/File: vbend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vbend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      FUNCTION VBEND (S)
C
C     VBEND  - compute bending energy levels
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C
C  Calls:
C     AITKF2 - Aitken interpolation
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     EBEND  - compute bending energy levels
C     LOCS   - locate position of s in grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(7)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /INTER1/ NINT
      LOGICAL LGS2(10)
      COMMON /LOGIC2/  LGS2
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2                       GCL0895
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      DATA IS /1/
C
      IF (S .LE. SS(1)) THEN
         FB = FBS(1)
         AB = QFBS(1)
         GB = GBS(1)
      ELSE IF (S .GE. SS(NSMAX)) THEN
         FB = FBS(NSMAX)
         AB = QFBS(NSMAX)
         GB = GBS(NSMAX)
      ELSE
         NP = NINT + 1
         CALL LOCS (IS, S)
         NN = NP/2 - 1
         IS = MAX(1, IS-NN)
         IS = MIN(IS, NSMAX-NP)
         FB = AITKF2(S, FBS(IS), F, SS(IS), NINT)
         AB = AITKF2(S, QFBS(IS), F, SS(IS), NINT)
         GB = AITKF2(S, GBS(IS), F, SS(IS), NINT)
      END IF
      IF (LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
         CALL COBEND (NBND1,NBDN2,FB,AB,GB,EB)
      ELSE
C  uncoupled bending energy level
         EB = EBEND (NBND1, FB, AB, GB)
         IF (NBND1 .EQ. NBND2) EB = 2.D0*EB
         IF (NBND1 .NE. NBND2) EB = EB + EBN (NBND2, FB, AB, GB)
      END IF
      VBEND = EB
C     WRITE(6, 6600) S, SS(1), SS(NSMAX), AA, VV
C6600 FORMAT(' S,SS(1),SS(NSMAX),AA,V=', 3F10.5, 1P4E15.7)
      RETURN
      END
