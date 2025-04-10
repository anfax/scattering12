!!!*************************************************************
! 文件/File: wkbtp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: wkbtp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WKBTP2 (E, SGN, X, XL, XR, IERR)
C
C     WKBTP2 - performs root search for turning points in potential for
C              u motion
C
C  Called by:
C     WKBTP  - find turning points in potential for u motion
C
C  Calls:
C     WKBPOT - potential for u motion for fixed s
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCONV
C
      IC = 0
C     WRITE (98, 9800) SGN, XL, XR, X
C9800 FORMAT (' WKBTP: SGN,XL,XR,X=', 1P4E13.5)
      LCONV = .FALSE.
   10 CONTINUE
         IC = IC + 1
         F = WKBPOT (X, DVDU, IRANGE) - E
         IF (IRANGE .NE. 0) THEN
            IF (SGN .GT. 0.0D0) THEN
               XR = MIN (XR, X)
            ELSE
               XL = MAX (XL, X)
            END IF
            X = 0.5D0*(XL + XR)
         ELSE
C              WRITE (98, 9801) IC, X, F, XL, XR
C9801          FORMAT (' IC,X,F,XL,XR=', I5, 1P4E15.7)
            LCONV = ABS(F) .LT. 1.D-8
            IF (.NOT. LCONV) THEN
               IF (SGN*F .GT. 0.0D0) THEN
                  XR = MIN (XR, X)
               ELSE
                  XL = MAX (XL, X)
               END IF
               IF (DVDU .NE. 0.0D0) THEN
                  X = X - F/DVDU
               ELSE
                  X = XL
               END IF
               IF (X .LE. XL .OR. X .GE. XR) X = 0.5D0*(XL + XR)
               LCONV = ABS(XR-XL) .LT. 1.D-9
            END IF
         END IF
      IF (.NOT.LCONV .AND. IC .LE. 100) GO TO 10
      IF (.NOT.LCONV) THEN
         WRITE (6, 601) E, X, F
         IERR = 3
      END IF
      RETURN
601   FORMAT ( ' WKBTP: ITERATION TO FIND TURNING POINT NOT CONVERGED,',
     *   ' E=', 1PE17.9, ', X=', E17.9, ', F=', E17.9)
      END
