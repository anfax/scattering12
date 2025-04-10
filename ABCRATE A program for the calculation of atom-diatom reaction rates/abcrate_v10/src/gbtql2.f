!!!*************************************************************
! 文件/File: gbtql2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: gbtql2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  GBTQL2
C***********************************************************************
C
      SUBROUTINE GBTQL2 (N,D,E,Z,IERR)
C
C     CALLED BY:
C                 GAUSSQ
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
C     METHOD, AND IS ADAPTED FROM THE EISPAK ROUTINE IMTQL2
C
C     ON INPUT=
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
C
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
C
C      ON OUTPUT=
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION MACHEP
      INTEGER I,J,K,L,M,N,II,MML,IERR
      DIMENSION D(N),E(N),Z(N)
C
C     ========== MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ==========
C
      MACHEP = 1.0D-16
C
      IERR = 0
      IF (N.EQ.1) GO TO 110
C
      E(N) = 0.0D0
      DO 70 L = 1, N
         J = 0
C
C     ========== LOOK FOR SMALL SUB-DIAGONAL ELEMENT ==========
C
   10    DO 20 M = L, N
            IF (M.EQ.N) GO TO 30
            IF (ABS(E(M)).LE.MACHEP*(ABS(D(M))+ABS(D(M+1)))) GO TO 30
   20    CONTINUE
C
   30    P = D(L)
         IF (M.EQ.L) GO TO 70
         IF (J.EQ.30) GO TO 100
         J = J+1
C
C     ========== FORM SHIFT ==========
C
         G = (D(L+1)-P)/(2.0D0*E(L))
         R = SQRT(G*G+1.0D0)
         G = D(M)-P+E(L)/(G+SIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M-L
C
C     ========== FOR I=M-1 STEP -1 UNTIL L DO -- ==========
C
         DO 60 II = 1, MML
            I = M-II
            F = S*E(I)
            B = C*E(I)
            IF (ABS(F).LT.ABS(G)) GO TO 40
            C = G/F
            R = SQRT(C*C+1.0D0)
            E(I+1) = F*R
            S = 1.0D0/R
            C = C*S
            GO TO 50
   40       S = F/G
            R = SQRT(S*S+1.0D0)
            E(I+1) = G*R
            C = 1.0D0/R
            S = S*C
   50       G = D(I+1)-P
            R = (D(I)-G)*S+2.0D0*C*B
            P = S*R
            D(I+1) = G+P
            G = C*R-B
C
C     ========== FORM FIRST COMPONENT OF VECTOR ==========
C
            F = Z(I+1)
            Z(I+1) = S*Z(I)+C*F
            Z(I) = C*Z(I)-S*F
C
   60    CONTINUE
C
         D(L) = D(L)-P
         E(L) = G
         E(M) = 0.0D0
         GO TO 10
   70 CONTINUE
C
C     ========== ORDER EIGENVALUES AND EIGENVECTORS ==========
C
      DO 90 II = 2, N
         I = II-1
         K = I
         P = D(I)
C
         DO 80 J = II, N
            IF (D(J).GE.P) GO TO 80
            K = J
            P = D(J)
   80    CONTINUE
C
         IF (K.EQ.I) GO TO 90
         D(K) = D(I)
         D(I) = P
C
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
C
   90 CONTINUE
C
      GO TO 110
C
C     ========== SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ==========
C
  100 IERR = L
  110 RETURN
C
C     ========== LAST CARD OF GBTQL2 ==========
C
      END
