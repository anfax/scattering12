!!!*************************************************************
! 文件/File: spl1b1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: spl1b1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SPL1B1(N,X,F,W,IJ,A,B,C,D)
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF THE INDEPENDENT VARIABLE
C           F = ORIGIN OF TABLE OF THE DEPENDENT VARIABLE
C           W = ORIGIN OF TABLE OF SECOND DERIVATIVES AS CALCULATED
C               BY SPL1D1
C          IJ = SPACING IN THE TABLES F AND W
C           A = ORIGIN OF TABLE OF THE CUBIC COEFFICIENT
C           B = ORIGIN OF TABLE OF THE QUADRATIC COEFFICIENT
C           C = ORIGIN OF TABLE OF THE LINEAR COEFFICIENT
C           D = ORIGIN OF TABLE OF THE CONSTANT COEFFICIENT
C
C
C     SPL1B1 CONVERTS THE SPLINE FIT DATA SUPPLIED BY SPL1D1
C     FROM THE FUNCTION AND ITS SECOND DERIVATIVE AT THE KNOTS,
C     TO THE FOUR COEFFICINTS OF THE CUBIC - A,B,C,D, WHERE
C      F IS NOW APPROXIMATED BY A*X**3 + B*X**2 + C*X + D
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(N),W(N),A(N),B(N),C(N),D(N)                      TCA0197
      SAVE C6, C3                                                       TCA1097
      DATA C6 /.1666 6666 6666D0/ , C3 /.3333 3333 3333D0/
C
C      GIVEN X,F,W COMPUTE A,B,C,D
C
      NM1 = N - 1
      DO 100 I = 1,NM1
      M = (I-1)*IJ + 1
      MPIJ = M + IJ
      IP1 = I + 1
      XIP1 = X(IP1)
      XI = X(I)
      H = XIP1 - XI
       FI = F(I)
      FIP1 = F(MPIJ)
      WIP1 = W(MPIJ)
      WI = W(M)
      AI = (WIP1 - WI)*C6
      T1 = XIP1*WI
      T2 = XI*WIP1
      BI = .5D0*(T1-T2)
      T1 = XIP1*T1
      T2 = XI*T2
      C(I) = (.5D0*(T2-T1) + FIP1 - FI)/H - AI*H
      A(I) = AI/H
      T1 = XIP1*T1
      T2 = XI*T2
      D(I) = (C6*(T1-T2) + FI*XIP1 - FIP1*XI)/H - C3*H*BI
      B(I) = BI/H
  100 CONTINUE
      RETURN
      E N D
