!!!*************************************************************
! 文件/File: wkbtp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: wkbtp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WKBTP (E, X1, X2, IERR)
C
C     WKBTP  - find turning points in potential for u motion
C
C  Called by:
C     PHSINT - compute phase integrals needed for WKB quantization
C
C  Calls:
C     WKBTP2 - performs root search for turning points in potential for
C              u motion
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /WKB2/   UMIN, VMIN, UMAX, VMAX
C
      IF (E-VMIN .LE. 1.D-8) THEN
         IERR = 1
      ELSE
         IF (E .GT. VMAX) THEN
            IERR = 2
            WRITE(6, 600) E, VMAX
         END IF
         IF (VMAX-E .LT. 1.D-8) THEN
            X2 = UMAX
         ELSE
            IERR = 0
            SGN = +1.0D0
            XL = UMIN
            XR = UMAX
            IF (X2 .LE. XL .OR. X2 .GE. XR) X2 = 0.5D0 * (XL + XR)
            CALL WKBTP2 (E, SGN, X2, XL, XR, IERR)
         END IF
         SGN = -1.0D0
         XL = -10000.0D0
         XR = UMIN
         IF (X1 .LE. XL .OR. X1 .GE. XR) X1 = 2.0D0*UMIN - X2
         CALL WKBTP2 (E, SGN, X1, XL, XR, IERR)
      END IF
      RETURN
  600 FORMAT ( ' WKBTP: E=', 1PE17.9, ', SHOULD BE LESS THAN VMAX=',
     *   E17.9)
      END
