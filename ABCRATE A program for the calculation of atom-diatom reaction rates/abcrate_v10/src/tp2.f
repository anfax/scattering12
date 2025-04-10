!!!*************************************************************
! 文件/File: tp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE TP2 (E, SR)
C
C     TP2    - find final turning point in adiabatic barrier in product
C              channel
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     CUBIC  - solve for roots of cubic polynomial
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LSET
      DIMENSION RRT(3), AIRT(3), SRT(4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
C     WRITE(6, 6600) E, NS2, NMID
C6600 FORMAT(' ENTERING TP2'/' E,NS2,NMID=', 1PE15.7, I5, I5)
C  check if the first point is higher than the energy
      IF (E .LT. VV2(NS2)) THEN
         S = SS(NS2)
         WRITE (6, 6000) S, E*CKCAL, VV2(NS2)*CKCAL
      ELSE
         IS = NS2
   10    CONTINUE
            IS = IS - 1
            LSET = E .LE. VV2(IS)
         IF (.NOT.LSET .AND. IS .GT. ISN) GO TO 10
         IF (LSET) THEN
C           WRITE(6, 6603) IS, E, SS(IS), SS(IS+1), VV2(IS),
C    *         VV2(IS+1)
C6603       FORMAT(' IS=', I5, ', E=', 1PE15.7, ', SS(IS),SS(IS+1)=',
C    *         0P2F15.7, ', VV2(IS),VV2(IS+1)=', 1P2E15.7)
C      S1 AND S2 ARE LEFT AND RIGHT BOUNDS ON THE TURNING POINT
            S1 = SS(IS)
            S2 = SS(IS+1)
            DD = DV2(IS) - E
C           WRITE(6, 6610) AV2(IS), BV2(IS), CV2(IS), DV2(IS), DD
C6610       FORMAT(' COEFFICIENTS'/1X, 1P5E15.7)
            CALL CUBIC (AV2(IS), BV2(IS), CV2(IS), DD, NREAL, RRT, AIRT)
            NRT = 0
            IF (NREAL .GE. 1) THEN
C              WRITE(6, 6604) NREAL, RRT
C6604          FORMAT(1X, I3, ' REAL ROOTS', 2X, 1P3E15.7)
               NRT = 0
               DO 20 I = 1,NREAL
                  IF (RRT(I) .LT. S1 .OR. RRT(I) .GT. S2) GO TO 20
                  NRT = NRT + 1
                  SRT(NRT) = RRT(I)
   20          CONTINUE
C              WRITE(6, 6605) NRT, SRT
C6605          FORMAT(1X, I3, ' ROOTS IN RANGE', 2X, 1P4E15.7)
            END IF
            IF (NRT .GT. 0) THEN
               S = SRT(NRT)
            ELSE
               S = 0.5D0 * (S1 + S2)
C              WRITE(6, 6606) S, S1, S2
C6606          FORMAT(' NO REAL ROOTS, S,S1,S2=', 3F15.7)
            END IF
         ELSE
            S = SS(ISN)
            WRITE (6, 6001) S, E*CKCAL, VV2(ISN)*CKCAL
         END IF
      END IF
      SR = S
      RETURN
C
6000  FORMAT(2X,T5,'Warning: In TP2 at the right most point s=', F10.5, 
     *             ', E=',1PE13.5, ', is less than V=', E13.5)
6001  FORMAT(2X,T5,'Warning: IN TP2 at the left most point s=', F10.5, 
     *             ', E=',1PE13.5, ', is greater than V=', E13.5)
      END
