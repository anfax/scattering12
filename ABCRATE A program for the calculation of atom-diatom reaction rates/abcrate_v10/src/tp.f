!!!*************************************************************
! 文件/File: tp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE TP (IFLG, E, SL, SR, SN, XT, LRFLC)
C
C     TP     - find turning points in adiabatic barrier
C
C  Called by:
C     LAG    - compute LAG probabilities
C     PSAG   - compute SAG-type probabilites
C
C  Calls:
C     Nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LRFLC, LRT
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
C
C     WRITE(6, 6600) IFLG, E, NS, LRT, NMID, LSYM
C6600 FORMAT(' ENTERING TP'/' IFLG,E,NS,LRT,NMID,LSYM=', I3, 1PE15.7,
C    *   I5, L2, I5, L2)
      SGN = 1.D0                                                        GCL1092
      ITP = 1
      IF (IFLG .EQ. 0) THEN
C  IFLG = 0, start from left most point
C  IFLG .NE. 0, start from last point
         IS0 = 1
         LRT = .FALSE.
C  check if the first point is higher than the energy, if so use the
C  first point as the turning point
         IF (E .LT. VV(1)) THEN
            SGN = -1.0D0
            WRITE (6, 6000) SV(1), E*CKCAL, VV(1)*CKCAL
            SL = SV(1)
            ITP = 2
            IFLG = 1
         END IF
      END IF
      IF (ITP .EQ. 1) THEN
         CALL TPSRCH (E, S, LRT, IFLG, SGN, LRFLC, IS0, ITP, SL)        GCL0896
         SL = S
         SGN = -1.0D0
      END IF
      IF (IFLG .NE. 0) THEN
         CALL TPSRCH (E, S, LRT, IFLG, SGN, LRFLC, IS0, ITP, SL)        GCL0896
         SR = S
         SN = (SL+SR)*0.5D0
         XT = (SR-SL)*0.5D0
      END IF
      RETURN
6000  FORMAT(2X,T5,'Warning: In TP the turning pt. at left most pt. ',
     *             's=',F10.5,', E=',1PE13.5,', is less than V=',E13.5)
      END
