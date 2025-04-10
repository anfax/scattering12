!!!*************************************************************
! 文件/File: cubic.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cubic.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE CUBIC(A,B,C,D,NREAL,RRT,AIRT)
C
C       SOLVE CUBIC EQUATION
C       A*X**3 + B*X**2 + C*X + D = 0
C
C       NREAL - NUMBER OF REAL ROOTS
C       RRT(1-3) - REAL PARTS OF THE 3 ROOTS
C       AIRT(1-3) - IMAG. PARTS OF THE 3 ROOTS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION RRT(3),AIRT(3)
      LOGICAL LFIRST
      SAVE LFIRST                                                       TCA1097
      DATA LFIRST/.TRUE./
      IF(LFIRST) THEN
         C3 = 1.D0/3.D0                                                 GCL1092
         RT3 = SQRT(3.D0)
         LFIRST = .FALSE.
      END  IF
      IF(A.EQ.0.D0) THEN
         IF(B.EQ.0.D0) THEN
            IF(C.EQ.0.D0) THEN
C   A,B,C = 0.,    NO ROOTS
               NREAL = 0
               DO 10 I = 1,3
                  RRT(I) = 0.D0                                         GCL1092
                  AIRT(I) = 0.D0                                        GCL1092
10             CONTINUE
            ELSE
C   A,B = 0,   REDUCES TO LINEAR EQUATION
               NREAL = 1
               RRT(1) = -D/C
               AIRT(1) = 0.D0                                           GCL1092
               DO 20 I = 2,3
                  RRT(I) = 0.D0                                         GCL1092
                  AIRT(I) = 0.D0                                        GCL1092
20             CONTINUE
            END IF
         ELSE
C   A = 0 ,  REDUCES TO QUADRATIC EQUATION
            AA = -0.5D0*C/B
            BB = D/B
            T2 = AA*AA - BB
            IF(T2.LT.0.D0) THEN                                         GCL1092
C      0 REAL ROOTS, 2 IMAG.
               NREAL = 0
               T = SQRT(-T2)
               RRT(1) = AA
               RRT(2) = AA
               RRT(3) = 0.D0                                            GCL1092
               AIRT(1) = T
               AIRT(2) = -T
               AIRT(3) = 0.D0                                           GCL1092
            ELSE
C      2 REAL ROOTS, 0 IMAG.
               NREAL = 2
               T = SQRT(T2)
               RRT(1) = AA + T
               RRT(2) = AA - T
               RRT(3) = 0.D0                                            GCL1092
               DO 30 I = 1,3
                  AIRT(I) = 0.0D0 
30             CONTINUE
            END IF
         END IF
      ELSE IF (D.EQ.0.0D0) THEN                                         GCL1092
         IF (C.EQ.0.0D0) THEN                                           GCL1092
            IF (B.EQ.0.0D0) THEN                                        GCL1092
C   3 real roots, all zero
               NREAL = 3
               DO 35 I = 1,3
                  RRT(I) = 0.0D0                                        GCL1092
                  AIRT(I) = 0.0D0                                       GCL1092
35             CONTINUE
            ELSE
C   3 real roots, two zero
               NREAL = 3
               DO 36 I = 1,2
                  RRT(I) = 0.0D0                                        GCL1092
                  AIRT(I) = 0.0D0                                       GCL1092
36             CONTINUE
               RRT(3) = -B/A
               AIRT(3) = 0.0D0                                          GCL1092
            END IF
         ELSE
C   one zero root plus roots of quadratic
            AA = -0.5D0*B/A
            BB = C/A
            T2 = AA*AA - BB
            IF(T2.LT.0.D0) THEN                                         GCL1092
C      1 REAL(ZERO) ROOTS, 2 IMAG.
               NREAL = 1
               T = SQRT(-T2)
               RRT(1) = 0.D0                                            GCL1092
               RRT(2) = AA
               RRT(3) = AA
               AIRT(1) = 0.D0                                           GCL1092
               AIRT(2) = T
               AIRT(3) = -T
            ELSE
C      3 REAL ROOTS
               NREAL = 3
               T = SQRT(T2)
               RRT(1) = 0.D0                                            GCL1092
               RRT(2) = AA + T
               RRT(3) = AA - T
               DO 38 I = 1,3
                  AIRT(I) = 0.0D0 
38             CONTINUE
            END IF
         END IF
      ELSE
C   Full cubic equation
         P = B/A
         Q = C/A
         R = D/A
         X0 = -C3*P
         P2 = P*P
         AA = Q - C3*P2
         BB = (2.D0*P2 - 9.D0*Q)*P/27.D0 + R
         IF(AA.EQ.0.D0) THEN
C     CASE 4 AA = 0, 1 REAL ROOT
            T = -BB
            X1 = SIGN(ABS(T)**C3,T)
            NREAL = 1
            RRT(1) = X1+X0
            AIRT(1) = 0.D0                                              GCL1092
            T = 0.5D0*X1
            RRT(2) = T+X0
            RRT(3) = T+X0
            AIRT(2) = RT3*X1
            AIRT(3) = -AIRT(2)
         ELSE IF(BB.EQ.0.D0) THEN
C   CASE 3,   BB=0
            RRT(1) = X0
            AIRT(1) = 0.D0                                              GCL1092
            IF(AA.LT.0.D0) THEN
C       FOR BB=0, AA<0,  3 REAL ROOTS
               NREAL = 3
               X1 = SQRT(-AA)
               RRT(2) = X1+X0
               RRT(3) = -X1+X0
               AIRT(2) = 0.D0                                           GCL1092
               AIRT(3) = 0.D0                                           GCL1092
            ELSE
C       FOR BB=0, AA>0,  1 REAL ROOT
               NREAL = 1
               RRT(2) = X0
               RRT(3) = X0
               X1 = SQRT(AA)
               AIRT(2) = X1
               AIRT(3) = -X1
            END IF
         ELSE
C
            BBH = 0.5D0*BB
            BBHS = BBH*BBH
            T2 = BBHS + (AA**3)/27.D0                                   GCL1092
            IF(T2.LT.0.D0) THEN
C   CASE 2,  T2< 0,    3 REAL ROOTS
               NREAL = 3
               T = SQRT(-T2)
               TT = -BBH
               THETA = C3*ATAN2(T,TT)
               R = SQRT(BBHS - T2)
               R3 = R**C3
               CS = COS(THETA)
               SN = SIN(THETA)
               X1 = R3*CS
               X2 = RT3*R3*SN
               RRT(1) = 2.D0*X1 + X0
               RRT(2) =-X1+X2+X0
               RRT(3) =-X1-X2+X0
               DO 40 I = 1,3
                  AIRT(I) = 0.D0                                        GCL1092
40             CONTINUE
            ELSE
C  CASE 1 T2 > 0. ,  ONE REAL ROOT
               NREAL = 1
               T = SQRT(T2)
               A3 = -BBH + T
               AP = SIGN(ABS(A3)**C3,A3)
               A3 = -BBH - T
               AM = SIGN(ABS(A3)**C3,A3)
               X1 = AM + AP
               X2 = -0.5D0*X1
               Y = 0.5D0*RT3*(AP-AM)
               RRT(1) = X1 + X0
               RRT(2) = X2 + X0
               RRT(3) = X2 + X0
               AIRT(1) = 0.D0                                           GCL1092
               AIRT(2) = Y
               AIRT(3) = -Y
            END IF
         END IF
      END IF
      IF(NREAL.GE.2) THEN
         NRM = NREAL - 1
         DO 170 I = 1,NRM
            IP = I + 1
            X1 = RRT(I)
            I1 = I
            DO 160 J = IP,NREAL
               IF(RRT(J).GT.X1) GO TO 160
               X1 = RRT(J)
               I1 = J
160         CONTINUE
            IF(I1.EQ.I) GO TO 170
            RRT(I1) = RRT(I)
            RRT(I) = X1
170      CONTINUE
      END IF
      RETURN
      END
