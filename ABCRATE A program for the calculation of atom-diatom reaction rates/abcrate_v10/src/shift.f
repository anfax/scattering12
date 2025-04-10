!!!*************************************************************
! 文件/File: shift.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: shift.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SHIFT (IS, ISHIFT, ISGN, SDIF, ISAD, NSSP)
C
C     SHIFT  - shift r.p. info in grid
C        ISHIFT>0, ISGN<0,  move IS-NSMAX to right ISHIFT locations
C        ISHIFT>0, ISGN>0,  move  1-IS    to right ISHIFT locations
C        ISHIFT<0, ISGN>0,  move  1-IS    to left  ISHIFT locations
C        ISHIFT<0, ISGN<0,  move IS-NSMAX to left  ISHIFT locations
C
C  Called by:
C     REFLEC - reflect MEP info across symmetry point
C     RPATH  - follow reaction path and store MEP info on grid
C     RPHSUM - summarize reaction path info
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NSSP(4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      IF (ISHIFT .LT. 0) THEN
         INC = 1
         IF (ISGN .LT. 0) THEN
            I2 = IS - 1
            I1 = I2 + ISHIFT
            N = NSDM - IS + 1
         ELSE
            I1 = 0
            I2 = -ISHIFT
            N = IS + ISHIFT
         END IF
      ELSE IF (ISHIFT .GT. 0) THEN
         INC = -1
         IF (ISGN .LT. 0) THEN
            I1 = NSDM + 1
            I2 = I1 - ISHIFT
            N = NSDM - IS + 1 - ISHIFT
         ELSE
            I2 = IS + 1
            I1 = I2 + ISHIFT
            N = IS
         END IF
      ELSE
         N = 0
      END IF
      IF (N .GT. 0) THEN
         DO 10 I = 1,N
            I1 = I1 + INC
            I2 = I2 + INC
            SS(I1) =SS(I2) + SDIF
            DS(I1) = DS(I2)
            VS(I1) = VS(I2)
            XKS(I1) = XKS(I2)
            X2(I1) = X2(I2)
            Y2(I1) = Y2(I2)
            UXS(I1) = UXS(I2)
            UYS(I1) = UYS(I2)
            FBS(I1) = FBS(I2)
            QFBS(I1) = QFBS(I2)
            GBS(I1) = GBS(I2)
            XMOMS(I1) = XMOMS(I2)
            CAPS(I1) = CAPS(I2)
   10    CONTINUE
         DO 20 I = 1,ISAD
            IF (ISGN*IS .GE. ISGN*NSSP(I)) NSSP(I) = NSSP(I) + ISHIFT
   20    CONTINUE
      END IF
      RETURN
      END
