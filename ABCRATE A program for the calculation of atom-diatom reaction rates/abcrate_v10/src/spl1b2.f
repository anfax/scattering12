!!!*************************************************************
! 文件/File: spl1b2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: spl1b2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SPL1B2(N,X,A,B,C,D,Y,TAB,IOP)
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF THE INDEPENDENT VARIABLE
C           A = ORIGIN OF TABLE OF CUBIC COEFFICIENTS
C           B = ORIGIN OF TABLE OF QUADRATIC COEFFICIENTS
C           C = ORIGIN OF TABLE OF LINEAR COEFFICIENTS
C           D = ORIGIN OF TABLE OF CONSTANT COEFFICIENTS
C           Y = THE POINT AT WHICH INTERPOLATION IS DESIRED
C         TAB = AN ARRAY OF DIMENSION 3 WHICH CONTAINS THE FUNCTION
C               VALUE, FIRST DERIVATIVE, AND SECOND DERIVATIVE AT Y
C         IOP = INTEGER SPECIFYING WHETHER DERIVATIVES ARE TO BE COMPUTED
C            .LE. 0 , F ONLY COMPUTED
C             = 1 , F AND FIRST DERIVATIVE ARE COMPUTED
C            .GE. 2 , F AND FIRST AND SECOND DERIVATIVES COMPUTED
C
C       THE FUNCTION F IS NOW APPROXIMATED BY THE CUBIC EQUATION
C          A*X**3 + B*X**2 + C*X + D
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),A(N),B(N),C(N),D(N),TAB(3)                         TCA0197
C
C     LOCATE Y IN THE X TABLE
C
      NM1 = N - 1
      IF(Y-X(1)) 10,10,20
   10 I = 1
      GO TO 30
   20 IF(Y-X(N)) 15,40,40
   40 I = NM1
      GO TO 30
   15 I = 1
      DO 25 K = 2,NM1
      IF(X(K).GT.Y) GO TO 30
   25 I = I + 1
   30 CONTINUE
C
C      CALCULATE F(Y)
C
      T = A(I)*Y
      TAB(1) = ((T + B(I))*Y + C(I))*Y + D(I)
C
C      CALCULATE THE FIRST DERIVATIVE OF F(Y)
C
      IF(IOP.LE.0) RETURN
      T = 3.D0*T
      TAB(2) = (T+2.D0*B(I))*Y + C(I)
C
C      CALCULATE THE SECOND DERIVATIVE OF F(Y)
C
      IF(IOP.EQ.1) RETURN
      TAB(3) = 2.D0*(T+B(I))
      RETURN
      E N D
