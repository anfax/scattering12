!!!*************************************************************
! 文件/File: d2dx2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: d2dx2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE D2DX2 (X0, Y0, F, ERR)
C
C     D2DX2 - compute matrix of second derivatives in Jacobi
C        coordinates.  Second derivative by taking numerical first
C        derivatives of analytic first derivatives.
C
C  Called by:
C     SADDLE - find saddle point(s) and do normal mode analysis
C     WELL   - print out info about local well in reaction path
C
C  Calls:
C     DERIV  - compute derivatives of potential w.r.t. x,y coordinates
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LCONV
      DIMENSION F(3) ,F0(4), F1(4), F2(4)
      EPSS0 = 1.D35
      EPS = 1.D-5
      DEL = 0.01D0
      DO 10 I = 1,4
         F2(I) = 0.0D0
   10 CONTINUE
      DO 40 I = 1,4
         DEL = DEL * 0.1D0
         X = X0 + DEL
         Y = Y0
         CALL DERIV(X, Y, DX, DY)
         X = X0 - DEL
         CALL DERIV(X, Y, D1, D2)
         F1(1) = 0.5D0*(DX-D1)/DEL
         F1(2) = 0.5D0*(DY-D2)/DEL
         Y = Y0 + DEL
         X = X0
         CALL DERIV(X, Y, DX, DY)
         Y = Y0-DEL
         CALL DERIV(X, Y, D1, D2)
         F1(3) = 0.5D0*(DY-D2)/DEL
         F1(4) = 0.5D0*(DX-D1)/DEL
         EPSS = 0.D0                                                    GCL1092
         DO 20 J = 1,4
            T = F1(J) - F2(J)
            F2(J) = F1(J)
            IF (F1(J) .NE. 0.0D0) T = T/F1(J)
            EPSS = EPSS + T*T
   20    CONTINUE
         IF (I .GT. 1) THEN
            LCONV = SQRT(EPSS).LT.EPS
            IF (LCONV) GO TO 50
            IF (EPSS .GT. EPSS0) GO TO 40
            EPSS0 = EPSS
         END IF
         DO 30 J = 1,4
            F0(J) = F1(J)
   30    CONTINUE
   40 CONTINUE
   50 CONTINUE
      IF (.NOT.LCONV) THEN
         DO 60 I = 1,4
            F1(I) = F0(I)
   60    CONTINUE
         EPSS = EPSS0
      END IF
      IF (ABS(F1(2)) .LT. 1.D-35 .AND. ABS((F1(2) - F1(4))/F1(2)) .GT.
     *   5.D0*EPS) WRITE (6, 601) EPS
      F(1) = F1(1)
      F(2) = 0.5D0*(F1(2) + F1(4))
      F(3) = F1(3)
      ERR = SQRT(EPSS)
      RETURN
601   FORMAT(/,2X,T5,'Warning: In D2DX2, F12 and F21 differ by more ',
     *               'than 5.0D0*',1PE9.3)
      END
