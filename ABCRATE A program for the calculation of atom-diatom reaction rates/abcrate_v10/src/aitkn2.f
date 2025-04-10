!!!*************************************************************
! 文件/File: aitkn2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: aitkn2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      SUBROUTINE AITKN2 (IS, S, X, Y, UX, UY)
C
C     AITKN2 - do Aitken interpolation of X, Y, UX, and UY from r.p.
C        grid
C
C  Called by:
C     LALFST - set up zeroth order path for LAG
C     LBETAS - solve for s of beta
C     LCSA   - compute cosine factors for LAG calculation
C
C  Calls:
C     AITKNF - Aitken interpolator
C     LOCS   - locate position of s in grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NFDIM=7)
      DIMENSION F(NFDIM)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /INTER1/ NINT
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL1096
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      IF (S .LE. SS(1)) THEN
C  If S is off the grid set values to last point on grid
         X = X2(1)
         Y = Y2(1)
         UX = UXS(1)
         UY = UYS(1)
      ELSE IF (S .GE. SS(NSMAX)) THEN
         X = X2(NSMAX)
         Y = Y2(NSMAX)
         UX = UXS(NSMAX)
         UY = UYS(NSMAX)
      ELSE
C  Interpolate
         NP = NINT + 1
         IF(NP.GT.NFDIM) THEN
            WRITE (6,6000) NINT,NFDIM-1
            STOP
         END IF
C  Locate grid point for S on the grid
         CALL LOCS (IS, S)
         NN = NP/2 - 1
         IS = MAX(1, IS-NN)
         IS = MIN(IS, NSMAX-NP+1)
C  Aitken interpolation of X,Y,UX,UY
         DO 10 I = 1,NP
            F(I) = X2(I+IS-1)
10       CONTINUE
         X = AITKNF(S, F, SS(IS), NINT)
         DO 20 I = 1,NP
            F(I) = Y2(I+IS-1)
20       CONTINUE
         Y = AITKNF(S, F, SS(IS), NINT)
         DO 30 I = 1,NP
            F(I) = UXS(I+IS-1)
30       CONTINUE
         UX = AITKNF(S, F, SS(IS), NINT)
         DO 40 I = 1,NP
            F(I) = UYS(I+IS-1)
40       CONTINUE
         UY = AITKNF(S, F, SS(IS), NINT)
         T = SQRT(UX*UX + UY*UY)
         UX = UX/T
         UY = UY/T
      END IF
      RETURN
6000  FORMAT(' ***** IN AITKN2, ORDER OF INTERPOLATION IS', I5,
     *   ' BUT MAXIMUM ALLOWED VALUE IS', I5)
      END
