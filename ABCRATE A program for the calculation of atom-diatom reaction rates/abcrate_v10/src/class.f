!!!*************************************************************
! 文件/File: class.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: class.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  CLASS
C***********************************************************************
C
      SUBROUTINE CLASS (KIND,N,ALPHA,BETA,B,A,MUZERO)
C
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
C        RECURRENCE RELATION
C
C             B P (X) = (X - A ) P   (X) - B   P   (X)
C              J J            J   J-1       J-1 J-2
C
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
C        AND THE ZERO-TH MOMENT
C
C             MUZERO = INTEGRAL W(X) DX
C
C        OF THE GIVEN POLYNOMIAL   WEIGHT FUNCTION W(X).  SINCE THE
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
C        GUARANTEED TO BE SYMMETRIC.
C
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
C        REQUIRE THE GAMMA FUNCTION.
C
C
C     CALLED BY:
C                 GAUSSQ
C     CALLS:
C                 DGAMA
C     ..................................................................
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION MUZERO
      DIMENSION A(N),B(N)
      SAVE PI                                                           TCA1097
      DATA PI / 3.141592653589793D0 /
C
      NM1 = N-1
      GO TO (10,30,50,70,90,110), KIND
C
C              KIND = 1=  LEGENDRE POLYNOMIALS P(X)
C              ON (-1, +1), W(X) = 1.
C
   10 MUZERO = 2.0D0
      DO 20 I = 1, NM1
         A(I) = 0.0D0
         ABI = I
         B(I) = ABI/SQRT(4*ABI*ABI-1.0D0)
   20 CONTINUE
      A(N) = 0.0D0
      RETURN
C
C              KIND = 2=  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
C
   30 MUZERO = PI
      DO 40 I = 1, NM1
         A(I) = 0.0D0
         B(I) = 0.5D0
   40 CONTINUE
      B(1) = SQRT(0.5D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 3=  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
C              ON (-1, +1), W(X) = SQRT(1 - X*X)
C
   50 MUZERO = PI/2.0D0
      DO 60 I = 1, NM1
         A(I) = 0.0D0
         B(I) = 0.5D0
   60 CONTINUE
      A(N) = 0.0D0
      RETURN
C
C              KIND = 4=  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
C              +INFINITY), W(X) = EXP(-X**2)
C
   70 MUZERO = SQRT(PI)
      DO 80 I = 1, NM1
         A(I) = 0.0D0
         B(I) = SQRT(I/2.0D0)
   80 CONTINUE
      A(N) = 0.0D0
      RETURN
C
C              KIND = 5=  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
C              BETA GREATER THAN -1
C
   90 AB = ALPHA+BETA
      ABI = 2.0D0+AB
      MUZERO = 2.0D0**(AB+1.0D0)*DGAMA(ALPHA+1.0D0)*DGAMA(BETA+1.0D0)/  GCL0992
     *   DGAMA(ABI)                                                     GCL0992
      A(1) = (BETA-ALPHA)/ABI
      B(1) = SQRT(4.0D0*(1.0D0+ALPHA)*(1.0D0+BETA)/((ABI+1.0D0)*ABI*ABI)
     *   )
      A2B2 = BETA*BETA-ALPHA*ALPHA
      DO 100 I = 2, NM1
         ABI = 2.0D0*I+AB
         A(I) = A2B2/((ABI-2.0D0)*ABI)
         B(I) = SQRT(4.0D0*I*(I+ALPHA)*(I+BETA)*(I+AB)/((ABI*ABI-1)*ABI*
     *      ABI))
  100 CONTINUE
      ABI = 2.0D0*N+AB
      A(N) = A2B2/((ABI-2.0D0)*ABI)
      RETURN
C
C              KIND = 6=  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
C              THAN -1.
C
  110 MUZERO = DGAMA(ALPHA+1.0D0)                                       GCL0992
      DO 120 I = 1, NM1
         A(I) = 2.0D0*I-1.0D0+ALPHA
         B(I) = SQRT(I*(I+ALPHA))
  120 CONTINUE
      A(N) = 2.0D0*N-1+ALPHA
      RETURN
      END
