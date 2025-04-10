!!!*************************************************************
! 文件/File: locs.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: locs.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LOCS (IS, S)
C
C     LOCS  - locate position of S in grid such that SS(IS) < or = S <
C        SS(IS+1)
C
C  Called by:
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR
C     AITKN2 - do Aitken interpolation of X, Y, UX, UY
C     DERS   - derivatives of morse turning point and zero pt. energy
C     VBEND  - compute bending energy levels
C     WKBSET - set up grid of WKB energy levels
C     VIBTA2 - compute vibrational period for u motion in product
C              channel
C     VSPLN2 - spline fit of adiabatic potential in product channel
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
C  check if first point is it
      IF (S .LT. SS(2)) THEN
         IS = 1
      ELSE IF (S .GE. SS(NSMAX)) THEN
         IS = NSMAX
      ELSE
C  search for IS, initial guess passed from calling routine
         IS = MAX(2, IS)
         IS = MIN(IS, NSMAX-1)
C  left and right bounds on IS
         ISL = 2
         ISR = NSMAX - 1
         IF (S .GE. SS(IS)) THEN
            ISL = IS
         ELSE
            ISR = IS
         END IF
C        IC = 0
C  search until left and right bounds differ by one
   10    CONTINUE
C  next guess assumes S on grid is linear with IS
C           IC = IC + 1
            IIS = IS + INT((S - SS(IS)) / (SS(IS) - SS(IS-1)))
            IIS = MAX(ISL, IIS)
            IIS = MIN(IIS, ISR)
            IF (IIS .EQ. IS) THEN
               IIS = (ISL + ISR)/2
            ELSE IF (IIS .EQ. ISL) THEN
               IIS = ISL + 1
            ELSE IF (IIS .EQ. ISR) THEN
               IIS = ISR - 1
            END IF
            IS = IIS
            IF (S .GE. SS(IS)) THEN
               ISL = IS
            ELSE
               ISR = IS
            END IF
C           IF (IC .GT. NSMAX) WRITE (6,6000) IC, IS, ISL, ISR, SS(IS),
C    *         S
C6000       FORMAT (' LOCS: IC, IS, ISL, ISR, SS(IS), S=', 4I5,
C    *      1P2E13.5)
C           IF (IC .GT. NSMAX+100) STOP
         IF (ISR-ISL .GT. 1) GO TO 10
         IS = ISL
      END IF
      RETURN
      END
