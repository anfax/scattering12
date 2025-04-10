!!!*************************************************************
! 文件/File: tpsrch.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tpsrch.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE TPSRCH (E, S, LRT, IFLG, SGN, LRFLC, IS0, ITP, SL)
C
C     TPSRCH - performs root search for turning points
C
C  Called by:
C     TP     - find turning points in adiabatic barrier
C
C  Calls:
C     CUBIC  - solve for roots of cubic polynomial
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LRFLC, LRT, LSET
      DIMENSION RRT(3), AIRT(3), SRT(4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
C  IFLG = 1,  successfully found tps
C  IFLG = 0,  no (more) tps found
      IFLG = 1
      IF (LRT) THEN
C        WRITE(6, 6601) SGN, NRC, NRT, SRT
C6601    FORMAT(' LRT=.T., SGN,NRC,NRT,SRT=', F3.1, 2I3, 4F15.7)
C  a root from the cubic polynomial left over, use it before searching
C     through the grid again
         NRC = NRC + 1
         S = SRT(NRC)
         IF(NRC.EQ.NRT) LRT = .FALSE.
      ELSE
         LRFLC = .FALSE.
C        WRITE(6, 6602) IS0, SGN
C6602    FORMAT(' LRT=.F., IS0=', I5, ', SGN=', F3.1)
         LSET = .FALSE.
C  search grid until E-V changes sign or end of grid is hit
         IS = IS0
   10    CONTINUE
            IF (LSET .OR. IS .GE. NS) GO TO 20
            IS = IS + 1
            LSET = (E-VV(IS))*SGN .LE. 0.D0                             GCL1092
         GO TO 10
   20    CONTINUE
         IS0 = IS
         IF (LSET) THEN
C  E-V changed sign between grid points IS-1 and IS
            IS = IS - 1
C           WRITE(6, 6603) IS, E, SV(IS), SV(IS+1), VV(IS), VV(IS+1)
C6603       FORMAT(' IS=', I5, ', E=', 1PE15.7, ', SV(IS),SV(IS+1)=',
C    *         0P2F15.7, ', VV(IS),VV(IS+1)=', 1P2E15.7)
C  S1 and S2 are left and right bounds on the turning point
            S1 = SV(IS)
            S2 = SV(IS+1)
            DD = DV(IS) - E
C           WRITE(6, 6610) AV(IS), BV(IS), CV(IS), DV(IS), DD
C6610       FORMAT(' COEFFICIENTS'/1X, 1P5E15.7)
C  V is expressed as a cubic spline, solve for zeros of cubic polynomial
            CALL CUBIC(AV(IS), BV(IS), CV(IS), DD, NREAL, RRT, AIRT)
            NRT = 0
            IF (NREAL .GE. 1) THEN
C              WRITE(6, 6604) NREAL, RRT
C6604          FORMAT(1X, I3, ' REAL ROOTS', 2X, 1P3E15.7)
C  readl roots found, check if they are between S1 and S2
               NRT = 0
               DO 30 I = 1,NREAL
                  IF (RRT(I) .LT. S1 .OR. RRT(I) .GT. S2) GO TO 30
                  NRT = NRT + 1
                  SRT(NRT) = RRT(I)
   30          CONTINUE
C              WRITE(6, 6605) NRT, SRT
C6605          FORMAT(1X, I3, ' ROOTS IN RANGE', 2X, 1P4E15.7)
            END IF
            IF (NRT .GT. 0) THEN
C  valid real roots found, use lowest one
               NRC = 1
               S = SRT(1)
               IF(NRT.GT.1) LRT = .TRUE.
            ELSE
C  no valid real root found, use the bound
               IF (ITP .EQ. 1) S = S1
               IF (ITP .EQ. 2) S = S2
C              WRITE(6, 6606) S, S1, S2
C6606          FORMAT(' NO REAL ROOTS, S,S1,S2=', 3F15.7)
            END IF
         ELSE
C  ran into end of grid
            IF (SGN .LT. 0.D0) THEN
C  looking for right turning point but couldn't find it
               IF (LSYM) THEN
C  for symmetric system just reflec the left turning point across
C     symmetry point
                  LRFLC = .TRUE.
                  S = -SL+2.0D0*SV(NS)
               ELSE
C  for nonsymmetric system write error message and use last point in
C     grid
                  WRITE(6, 6001) SV(NS), E*CKCAL, VV(NS)*CKCAL
                  S = SV(NS)
               END IF
            ELSE
C  looding for left turning point but couldn't find it, return with
C     IFLG=0 (signals end of tps)
               IFLG = 0
            END IF
         END IF
      END IF
      RETURN
C
6001  FORMAT(2X,T5,'Warning: In TPSRCH the turning pt. at right most ',
     *     'point s=',F10.5,', E=',1PE13.5, ', is less than V=', E13.5)
      END
