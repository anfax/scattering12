!!!*************************************************************
! 文件/File: liadib.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: liadib.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LIADIB (S,E,ALF,BETA,X0,Y0,UX,UY,XINT,VAD,COSCHI,
     *   LNONAD)
C
C     LIADIB -  evaluate adiabatic integrand
C
C  Called by:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C
C  Calls:
C     LALF0  - computes x and y as function of beta
C     SPL1B2 - evaluate spline fit
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LNONAD
      DIMENSION TAB(3)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /ADIAB1/ VADX(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /THLAG1/ XM1, YM1, X1, Y1, UXM1, UYM1, UX1, UY1, SL, SR,
     *XDIF, XSUM, YDIF, YSUM
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
C    compute contribution to phase integral from VA region
      S1 = S
      IF (.NOT.LNONAD .AND. LSYM .AND. S .GE. SS(NMID)) THEN
         S1 = 2.D0 * SS(NMID) - S                                       GCL1092
      END IF
      IF (LNONAD) THEN
C  'product-state' adiabatic potential
         CALL SPL1B2(NS2, SS, AV2, BV2, CV2, DV2, S1, TAB, 1)
      ELSE
C  'reactant-state' adiabatic potential
         CALL SPL1B2(NS, SV, AV, BV, CV, DV, S1, TAB, 1)
      END IF
      XINT = 0.0D0                                                      GCL1092
      COSCHI = 0.0D0                                                    GCL1092
      VAD = TAB(1)
      T = VAD - E
      IF (T .GT. 0.D0) THEN                                             GCL1092
         XINT = SQRT(T)
C    compute cos(chi),  first deriv. of xab,yab w.r.t. beta
C    compute XAB, YAB for beta +/- 1.D-4
         STEPB = 1.D-4                                                  GCL1092
         BET = BETA - STEPB
         CALL LALF0 (BET,X0,Y0)
         XABT1 = (1.D0 - ALF) * X0 + ALF * (XSUM + BET * XDIF)          GCL1092
         YABT1 = (1.D0 - ALF) * Y0 + ALF * (YSUM + BET * YDIF)          GCL1092
         BET = BETA + STEPB
         CALL LALF0 (BET,X0,Y0)
         XABT2 = (1.D0 - ALF) * X0 + ALF * (XSUM + BET * XDIF)          GCL1092
         YABT2 = (1.D0 - ALF) * Y0 + ALF * (YSUM + BET * YDIF)          GCL1092
         DX = (XABT2 - XABT1) * 0.5D0 / STEPB                           GCL1092
         DY = (YABT2 - YABT1) * 0.5D0 / STEPB                           GCL1092
         COSCHI =  ABS(DX*UY - DY*UX)
         XINT = XINT*COSCHI
      END IF
      RETURN
      END
