!!!*************************************************************
! 文件/File: grad.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: grad.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE GRAD (XP, YP, DX, DY, STEP, DELS, DELSS, DELSV, VNOW,
     *   NSTAB, IRET)
C
C     GRAD   - follow gradient
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     VMIN   - find minimum energy along u coordinate
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSET
C  follow gradient
      IRET = 1
      ERREPS = 0.01D0*STEP
      DELSV2 = DELSV
C  loop over segments of NSTAB steps until the total distance is DELSV
C  (IRET=1) or AC<0.9 or NSTAB=1 (IRET=2)
   10 CONTINUE
         X = XP
         Y = YP
         XSV = X
         YSV = Y
         DXSV = DX
         DYSV = DY
         DELS = 0.D0                                                    GCL1092
C  loop over NSTAB steps of length STEP
         I = 0
   20    CONTINUE
            I = I + 1
            X0 = X
            Y0 = Y
            DX0 = DX
            DY0 = DY
            X = X - STEP*DX
            Y = Y - STEP*DY
            DELS = DELS + STEP
            DELSS = DELSS + STEP
            CALL DERIV(X, Y, DX, DY)
            DMAG = SQRT(DX*DX + DY*DY)
            DX = DX/DMAG
            DY = DY/DMAG
C   Check for oscillations of the gradient.  AC is the cos of the angle
C   between the old and new gradient vectors.  Restricting it to be
C   greater than .9 is equivalent to restricting the angle to 25.8
C   degrees.
            AC = DX*DX0 + DY*DY0
            IF (AC .LT. 0.9D0) IRET = 2
         IF (I .LT. NSTAB .AND. DELSS .LT. DELSV2 .AND. IRET .NE. 2)
     *      GO TO 20
         IF (IRET .EQ. 2)  GO TO 30
C  NSTAB steps successfully completed, now stabilize
         XP = X
         YP = Y
         DXP = DX
         DYP = DY
         CALL VMIN(XP, YP, DXP, DYP, STEP, VNOW, .TRUE., IERR)
         LSET = .TRUE.
         IF (IERR .NE. 0) THEN
C  stabilization unsuccessfull, try again with fixed gradient vector
            XP = X
            YP = Y
            DXP = DX
            DYP = DY
            CALL VMIN(XP, YP, DXP, DYP, STEP, VNOW, .FALSE., IERR)
            IF (IERR .NE. 0) THEN
C  still couldn't stabilize, write error message and forge ahead
               WRITE (6, 6000) NSTAB, DELSS, X, Y, DX, DY
               XP = X
               YP = Y
C  DELSV2 is reset so that gradient following will continue until out
C     of the region where trouble is
               IF (DELSS .GE. DELSV2) DELSV2 = DELSV2 + STEP*NSTAB
               LSET = .FALSE.
            END IF
         END IF
         IF (LSET) THEN
            DX = DXP
            DY = DYP
            T1 = XP - X0
            T2 = YP - Y0
            DELSS = DELSS + SQRT(T1*T1+T2*T2) - STEP
C  error estimate
            T1 = XP - X
            T2 = YP - Y
            ERR = SQRT(T1*T1+T2*T2)
            IF (ERR .GT. STEP) NSTAB = MAX(1, NSTAB-10)
            IF (NSTAB .EQ. 1) THEN
               IRET = 2
            ELSE
               IF (ERR .LT. ERREPS) NSTAB = NSTAB + 5
            END IF
         END IF
      IF (DELSS .LT. DELSV2 .AND. IRET .NE. 2) GO TO 10
   30 CONTINUE
      IF (IRET .EQ. 2) THEN
C  oscillations were too great, reset to saved values and return
         XP = XSV
         YP = YSV
         DX = DXSV
         DY = DYSV
         DELSS = DELSS - DELS
         DELS = 0.D0                                                    GCL1092
      END IF
      RETURN
6000  FORMAT(1X,T5,'Warning: In GRAD there is a problem with VMIN; ',
     *             'NSTAB = ',I10,' , DELSS = ',1PE13.5,
     *       /,1X,T14,'X, Y, DX, DY = ',4(1PE13.5,1X))
      END
