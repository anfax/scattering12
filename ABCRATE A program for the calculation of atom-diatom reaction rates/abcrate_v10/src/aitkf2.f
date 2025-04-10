!!!*************************************************************
! 文件/File: aitkf2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: aitkf2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      FUNCTION AITKF2 (Y, F, FX, X, NINT)
C
C     AITKF2 - checks for "bad" values before calling Aitken
C        interpolator
C
C  Called by:
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR
C     KAPVA  - compute kappas
C     PLAG   - compute primitive LAG probabilities from thetas
C     VBEND  - compute bending energy levels
C     WKBSET - set up grid of WKB energy levels
C
C  Calls:
C     AITKNF - Aitken interpolator
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FX(*), F(*), X(*)                                       GCL1096
C
      IF (NINT .LT. 0) THEN
         AITFK2 = 0.0D0
      ELSE IF (NINT .EQ. 0) THEN
         AITFK2 = F(1)
      ELSE
         N = NINT
         I0 = 0
         TOLD = 0.D0                                                    GCL1092
C  Loop over order of interpolation, N.  Start at N=NINT and decrease
C     by one until successful interpolation or N=1
         ICOUNT = 0
   10    CONTINUE
            ICOUNT = ICOUNT + 1
            IF (ICOUNT .GT. 100) THEN
               WRITE (6, 6001) NINT, N
               WRITE (6, 6000) I0, N, NP1, Y, (X(I0+I), F(I0+I),
     *            I=1,NP1)
               STOP 'AITKF2 1'
            END IF
            NP1 = N + 1
C  Loop over points and check that no X values are identical and set
C     scratch array FX = F
            DO 20 I = 1,NP1
               IF (I .GT. 1) THEN
                  IF (X(I+I0). EQ. X(I+I0-1)) THEN
                     WRITE (6, 6000) I0, N, NP1, Y, (X(I0+J), F(I0+J),
     *                  J=1,NP1)
                      STOP 'AITKF2 2'
                   END IF
               END IF
               FX(I) = F(I+I0)
   20       CONTINUE
C  Aitken interpolation
            T = AITKNF(Y, FX, X(I0+1), N)
            IF (N. LE. 1) GO TO 50
C  Check that interpolated value is in correct range
            FMIN = F(I0+1)
            FMAX = F(I0+1)
            DO 30 I = 2,NP1
               FMIN = MIN(FMIN, F(I+I0))
               FMAX = MAX(FMAX, F(I+I0))
   30       CONTINUE
            IF (T. GE. FMIN .AND. T .LE. FMAX) GO TO 50
C  Interpolated value out of range; check that the difference in two
C     successive interpolation is less that 0.1
            TT = TOLD - T
            IF (T. NE. 0.D0) TT = TT/T
            IF (ABS(TT) .LT. 0.1D0) GO TO 50
C  Interpolation not successful; decrease N and try again.
            TOLD = T
            N = N - 1
            NH = NP1/2
            IF (X(NH+I0+1) .LT. Y) I0 = I0 + 1
         GO TO 10
   50    CONTINUE
         AITKF2 = T
      END IF
      RETURN
6000  FORMAT(/,2X,T5,'In AITKF2: I0, N, NP1, and Y = ',3(I5,1X),1PE15.7,
     *       /,2X,T16,'X, F = ',2(1PE15.7,1X))
6001  FORMAT(/,2X,T5,'Error in AITKF2: NINT and N = ',2(I5,1X))
      END
