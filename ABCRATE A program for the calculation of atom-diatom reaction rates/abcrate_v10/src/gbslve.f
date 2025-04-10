!!!*************************************************************
! 文件/File: gbslve.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: gbslve.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  GBSLVE
C***********************************************************************
C
      FUNCTION GBSLVE (SHIFT,N,A,B)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
C
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,
C
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
C       THE N-TH POSITION.
C
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
C
C       CALLED BY:
C                  GAUSSQ
C                                
      DIMENSION A(N),B(N)
C
      ALPHA = A(1)-SHIFT
      NM1 = N-1
      DO 10 I = 2, NM1
         ALPHA = A(I)-SHIFT-B(I-1)**2/ALPHA
   10 CONTINUE
      GBSLVE = 1.0D0/ALPHA
      RETURN
      END
