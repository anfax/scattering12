!!!*************************************************************
! 文件/File: root2d.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: root2d.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE ROOT2D (FUNC, XX, YY, DEL, EPS1, EPS2, NMAX)
C
C     ROOT2D - 2d root search
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C     SADDLE - find saddle point(s) and do normal mode analysis
C
C  Calls:
C     FUNC   - function passed as agrument to routine
C
C     Solve F(X,Y) = G(X,Y) = 0.
C     Method:
C          F(X,Y) and G(X,Y) are fit to linear equations in X,Y
C            F = A1 + B1*X + C1*Y
C            G = A2 + B2*X + C2*Y
C          The coefficients are computed by evaluating F and G
C          at three different points. The resulting three
C          simultaneous linear equations can then be solved.
C          Given the coefficents the equations are extrapolated
C          to zero for a new guess of the root (X,Y).
C          This procedure is iterated to convergence.
C
C     Paramters:
C          FUNC  -  user supplied subroutine FUNC(X, Y, F, G) which
C                   evaluates F and G for given values of X,Y.
C          XX,YY  - upon input these are the initial guess for the
C                   root.  Upon output they contain the root.
C          DEL -    step increment to generate the necessary second
C                   and third points for initialization of the
C                   procedure.  The three initial points are
C                   (X,Y),(X+DEL,Y),(X,Y+DEL)
C          EPS1  -  convergence criteria number one.  If the absolute
C                   values of both F and G are less than EPS1, control
C                   is returned to the calling routine.
C          EPS2  -  convergence criteria number two.  If both X and Y
C                   change by less than EPS2 on an iteration, control is
C                   returned to the calling routine.
C          NMAX  -  if number of iterations exceed NMAX an error is
C                   printed, NMAX is set to zero, and control is
C                   returned to the calling routine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FUNC                                                     GCL1096
      DIMENSION X(3), Y(3), F(3), G(3)
C   machine dependent constant, largest real value
      SAVE RMAX                                                         TCA1097
      DATA RMAX /1.D38/
      N = 0
      X(1) = XX
      Y(1) = YY
      CALL FUNC(X(1), Y(1), F(1), G(1))
      IF (ABS(F(1)) .LT. EPS1 .AND. ABS(G(1)) .LT. EPS1) RETURN
      X(2) = XX+DEL
      Y(2) = YY
      CALL FUNC(X(2), Y(2), F(2), G(2))
      X(3) = XX
      Y(3) = YY + DEL
   10 CONTINUE
         CALL FUNC(X(3), Y(3), F(3), G(3))
         IF (ABS(F(3)) .LT. EPS1 .AND. ABS(G(3)) .LT. EPS1) GO TO 40
         SUM1 = 0.D0                                                    GCL1092
         SUM2 = 0.D0                                                    GCL1092
         SUM3 = 0.D0                                                    GCL1092
         DO 20 I = 1,3
            IP = MOD(I, 3) + 1
            IPP  = MOD(IP, 3) + 1
            T = F(IP)*G(IPP) - G(IP)*F(IPP)
            SUM1 = SUM1 + T
            SUM2 = SUM2+ T*X(I)
            SUM3 = SUM3 + T*Y(I)
   20    CONTINUE
         DO 30 I = 1,2
            IP = I + 1
            X(I) = X(IP)
            Y(I)= Y(IP)
            F(I) = F(IP)
            G(I) = G(IP)
   30    CONTINUE
         IF (SUM1 .EQ. 0.0D0) THEN
            WRITE(6, 601) N, X(3), Y(3)
            NMAX=0
            GO TO 40
         END IF
         X(3) = SUM2/SUM1
         Y(3) = SUM3/SUM1
         IF (ABS(X(3)-X(2)) .LT. EPS2 .AND. ABS(Y(3)-Y(2)) .LT. EPS2)
     *      GO TO 40
         N = N + 1
         IF (N .GE. NMAX) THEN
            WRITE(6, 600) NMAX, X(3), Y(3)
            NMAX = 0
            GO TO 40
         END IF
      GO TO 10
   40 CONTINUE
      IF (X(3) .LT. RMAX .AND. Y(3) .LT. RMAX) THEN
         XX = X(3)
         YY = Y(3)
      ELSE
         WRITE(6, 603) N, X(3), Y(3)
         NMAX = 0
      END IF
      RETURN
  600 FORMAT (1X,T5,'Warning: In ROOT2D 2D search not converged, ',
     *              'NMAX = ', I5, 
     *     /,1X,T14,'last X, Y = ', 1PE13.5, ', ', 1PE13.5) 
  601 FORMAT (1X,T5,'Warning: In ROOT2D - division by zero, ITER = ',I5, 
     *     /,1X,T14,'last X, Y = ', 1PE13.5, ', ', 1PE13.5)  
  603 FORMAT (1X,T5,'Warning: In ROOT2D - overflow, ITER = ', I5, 
     *     /,1X,T14,'last X, Y = ', 1PE13.5, ', ', 1PE13.5)
      END
