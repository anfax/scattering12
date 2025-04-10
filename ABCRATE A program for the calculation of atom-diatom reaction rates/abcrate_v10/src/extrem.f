!!!*************************************************************
! 文件/File: extrem.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: extrem.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE EXTREM (IS, FOLD, FNEW, LIND, SGN, NMAX, NEXTR, ISEXTR,
     *   FEXTR)
C
C     EXTREM - finds extrema on grid and stores locations
C
C  Called by:
C     GTST   - compute free energies, CVT and ICVT rates
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LIND
      DIMENSION ISEXTR(4), FEXTR(4)
C
C  LIND should be false for first point on grid
      IF (.NOT.LIND) THEN
         SGN = FNEW - FOLD
         IF (SGN .NE. 0.D0) THEN
            IF (SGN .LT. 0.D0) THEN
C  asymptotic reactant value is a max
               NEXTR = 1
               FEXTR(1) = FOLD
               ISEXTR(1) = 0
            ELSE
               NEXTR = 0
            END IF
            LIND = .TRUE.
            SGN = SIGN(1.0D0, SGN)
         END IF
      ELSE
         T1 = FNEW-FOLD
         IF (T1*SGN .LT. 0.0D0) THEN
            SGN = -SGN
            NEX = NEXTR+1
            NEXTR = NEX
            IF (SGN .GT. 0.0D0) THEN
C  extremum was a min
               IF (NEX .LE. 4) THEN
                  FEXTR(NEX) = FOLD
                  ISEXTR(NEX) = IS - 1
               ELSE
                  IF (FEXTR(4) .GE. FOLD) THEN
                     FEXTR(4) = FOLD
                     ISEXTR(4) = IS-1
                  END IF
               END IF
            ELSE
C  extremum was a max
               IF (NEX .LE. 3) THEN
                  FEXTR(NEX) = FOLD
                  ISEXTR(NEX) = IS-1
                  NMAX = 1
                  IF (FEXTR(1) .GT. FOLD) NMAX = 3
               ELSE
                  IMAX = NMAX
                  IF (FEXTR(IMAX) .LE. FOLD) THEN
                     IF (IMAX .EQ. 1) THEN
                        FEXTR(1) = FEXTR(3)
                        ISEXTR(1) = ISEXTR(3)
                        FEXTR(2) = FEXTR(4)
                        ISEXTR(2) = ISEXTR(4)
                     ELSE
                        IF (FEXTR(2) .LE. FEXTR(4)) THEN
                           FEXTR(2) = FEXTR(4)
                           ISEXTR(2) = ISEXTR(4)
                        END IF
                     END IF
                     FEXTR(3) = FOLD
                     ISEXTR(3) = IS-1
                     FEXTR(4) = +1.D35
                     NMAX = 1
                     IF (FEXTR(1) .GT. FOLD) NMAX = 3
                  END IF
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
