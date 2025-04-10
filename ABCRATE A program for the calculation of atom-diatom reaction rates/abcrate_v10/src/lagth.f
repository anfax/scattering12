!!!*************************************************************
! 文件/File: lagth.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lagth.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LAGTH (E,ALF,THETA,IC3D,PERC,IERR,LNONAD,NFINAL)
C
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C
C  Called by:
C     LALFMN - find minimum theta(alf)
C
C  Calls:
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR
C     FUMCP  - compute outer turning point for u motion (Marcus-Coltrin
C              path)
C     LALF0  - computes x and y as function of beta
C     LBETAS - solve for s of beta
C     LIADIB -  evaluate adiabatic integrand
C     PEF    - evaluate potential
C     VBEND  - compute bending energy levels
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LNONAD,LEXAD
      DIMENSION THETA(2), PERC(2), SUM(2)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      DIMENSION BETAPT(NQKPDM),WEIGHT(NQKPDM,2),XINTG(NQKPDM)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /THLAG1/ XM1, YM1, X1, Y1, UXM1, UYM1, UX1, UY1, SL, SR,
     *XDIF, XSUM, YDIF, YSUM
C
      IERR = 0
      DELSMX = DELSV
      DELSMN = 0.001D0*DELSV                                            GCL1092
C zeta is the progress variable along the LCG path used in the POLYRATE
C   code and goes from 0 to zetamx defined below.  beta is the progress
C   variable used here and goes from -1 to 1.
C      IDBG = 97
C      ZETAMX = SQRT(XDIF*XDIF+YDIF*YDIF)
C      WRITE (IDBG, 9700) E*CKCAL,ALF,XM1,YM1,X1,Y1,ZETAMX
C 9700 FORMAT(1X,131(1H-)/' LAGTH, E,ALF,=', 1P2E15.7,/,
C     *   ' XM1,YM1,X1,Y1,ZETAMX=',1P5E15.7)
      DO 10 I = 1,2
         THETA(I) = 0.D0                                                GCL1092
         PERC(I) = 0.D0                                                 GCL1092
         SUM(I) = 0.D0                                                  GCL1092
10    CONTINUE
      IB = 0
      DBET = 1.0D0/DBLE(NSEG3)                                          GCL1092
      BSUM = -1.0D0 - DBET                                              GCL1092
      DO 15 I = 1,NSEG3
         BSUM = BSUM + 2.0D0*DBET                                       GCL1092
         DO 14 IQ = 1,NQ32
            IB = IB + 1
            BETAPT(IB) = BSUM + DBET*PT3(IQ)
            WEIGHT(IB,1) = WT3(IQ,1)*DBET
            WEIGHT(IB,2) = WT3(IQ,2)*DBET
14       CONTINUE
15    CONTINUE
      NBETA = NSEG3*NQ32
C ----------------------------------------------------------------------
C Loop over two directions of integration - from reactant (IDIR=1) and
C   product (IDIR=-1) sides.
      IDIR = 1
20    CONTINUE
         IF (IDIR.EQ.1) THEN
            XAB = XM1
            YAB = YM1
            UX = UXM1
            UY = UYM1
            S = SL
            ISINT = 1
            IB = 0
            LEXAD = .FALSE.
C            WRITE (IDBG, 9701)
C9701        FORMAT (' REACTANT DIRECTION, ADIABATIC REGION')
         ELSE
C   product side
            XAB = X1
            YAB = Y1
            UX = UX1
            UY = UY1
            S = SR
            ISINT = NSMAX
            IB = NBETA + 1
            LEXAD = LNONAD
C            WRITE (IDBG, 9702)
C9702        FORMAT (' PRODUCT DIRECTION, ADIABATIC REGION')
         END IF
C         WRITE (IDBG,9703)
C9703     FORMAT(4X,'IB',2X,'BETA',6X,'ZETA',7X,'S',9X,'VAD',7X,'COS',
C     *      6X,'XINT',7X,'X',9X,'Y',9X,'X(S)',6X,'Y(S)',6X,'UALF',
C     *      6X,'UMCP')
C.......................................................................
C Loop over grid points in adiabatic region
30       CONTINUE
            IB = IB + IDIR
            BETA = BETAPT(IB)
C   solve for s of beta
            CALL LBETAS (E,ALF,BETA,XAB,YAB,X,Y,UX,UY,UALF,
     *       DELSMN,DELSMX,S,ISINT,IERR)
            IF (IERR.NE.0) THEN
               WRITE (6,6000) ALF,BETA
               STOP 'LAGTH 1'
            END IF
C    compute UMCP
            CALL AITKEN (ISINT, S, VXX, D, XK, CP)
            SGN = 1.0D0                                                 GCL1092
            UTP = FUMCP (S, D, XK, SGN, LEXAD,NFINAL)
            UKAP = UTP
            IF (CP.GT.0.0D0) UKAP = MAX(0.0D0, 1.0D0/CP)                GCL1092
            UMCP = MIN(UTP,UKAP)
C    Check if inside adiabatic region (UALF<=UMCP)
            IF (UALF.LE.UMCP) THEN
C   evaluate adiabatic integrand
               CALL LIADIB (S,E,ALF,BETA,X0,Y0,UX,UY,XINT,VAD,COSCHI,
     *            LEXAD)
C               WRITE (IDBG,9704) IB,BETA,(1.0D0+BETA)*ZETAMX,          GCL1092
C     *            S,VAD,COSCHI/ZETAMX,XINT,XAB,YAB,X,Y,UALF,UMCP
C9704           FORMAT(1X,I5,1P12E10.2)
               IF(IDIR.EQ.1) THEN
C First pass from reactant side, save integrand in array
                  XINTG(IB) = XINT
               ELSE 
C Second pass, check if grid point has already been computed from
C    reactant side and if so, select minimum
                  IF (IB.LE.IB0-1) XINT = MIN(XINT,XINTG(IB))
                  XINTG(IB) = XINT
               END IF
C
               IF ((IDIR .GT. 0 .AND. IB .LT. NBETA) .OR. 
     *             (IDIR .LT. 0 .AND. IB .GT. 1)) GO TO 30
            END IF
C.......................................................................
C Integration point outside adiabatic region or at end of integration
C    grid
         IF(IDIR.EQ.1)THEN
            IDIR = -1
            XAB0 = XAB
            YAB0 = YAB
            SB0 = S
            IF (UALF.LE.UMCP) THEN
C Entire integration range in the adiabatic region
               IB0 = NBETA+1
               BETA0 = 1.0D0                                            GCL1092
            ELSE
               IB0 = IB
               BETA0 = BETA
C               WRITE (IDBG,9710) IB,BETA,(1.0D0+BETA)*ZETAMX,          GCL1092
C     *            S,XAB,YAB,X,Y,UALF,UMCP
C9710           FORMAT(1X,I5,1P3E10.2,30X,1P6E10.2)
            ENDIF
            GO TO 20
         END IF
C ----------------------------------------------------------------------
      XAB1 = XAB
      YAB1 = YAB
      SB1 = S
      IF (UALF.LE.UMCP) THEN
         IB1 = 0
         BETA1 = -1.0D0                                                 GCL1092
      ELSE 
         IB1 = IB
         BETA1 = BETA
C         WRITE (IDBG,9710) IB,BETA,(1.0D0+BETA)*ZETAMX,                GCL1092
C     *      S,XAB,YAB,X,Y,UALF,UMCP
      END IF
C
C BETA0 and BETA1 are the limits on the nonadiabatic region and IB0 and
C    IB1 are the indices of the first grid points in the nonadiabatic
C    region.
      IF (IB1 .GE. IB0) THEN
C Sum up contributions from adiabatic region
         IF (IB0.GT.1) THEN
            DO 40 IB = 1,IB0-1
               SUM(1) = SUM(1) + WEIGHT(IB,1)*XINTG(IB)
               SUM(2) = SUM(2) + WEIGHT(IB,2)*XINTG(IB)
40          CONTINUE
         END IF
         IF (IB1.LT.NBETA) THEN
            DO 50 IB = IB1+1,NBETA
               SUM(1) = SUM(1) + WEIGHT(IB,1)*XINTG(IB)
               SUM(2) = SUM(2) + WEIGHT(IB,2)*XINTG(IB)
50          CONTINUE
         END IF
C
C Do the integral in the nonadiabatic region.
C
C Compute the correction term for the effective potential in the
C    nonadiabatic region at the adiabatic/nonadiabatic boundary
         DB = BETA1 - BETA0
         IF (DB.GT.0.0D0) THEN
            DXDB = (XAB1 - XAB0)/DB
            DYDB = (YAB1 - YAB0)/DB
            DXIDB = SQRT(DXDB*DXDB + DYDB*DYDB)
         ELSE
            STEPB = 1.D-4                                               GCL1092
            BET = BETA0 - STEPB
            CALL LALF0 (BET,X0,Y0)
            XABT1 = (1.D0 - ALF) * X0 + ALF * (XSUM + BET * XDIF)       GCL1092
            YABT1 = (1.D0 - ALF) * Y0 + ALF * (YSUM + BET * YDIF)       GCL1092
            DXDB = (XABT1 - XAB0)/STEPB
            DYDB = (YABT1 - YAB0)/STEPB
            DXIDB = SQRT(DXDB*DXDB + DYDB*DYDB)
         END IF
         IF (IC3D .EQ. 2) THEN
C   compute bending contribution
            VB0 = VBEND(SB0)
            DVB = VBEND(SB1) - VB0
         END IF
C         WRITE (IDBG, 9707) VB0,VB0+DVB
C9707     FORMAT(' MIDDLE POTENTIAL REGION, VCORR(S0),VCORR(S1)=',
C     *      1P2E15.7,/,4X,'IB',2X,'BETA',6X,'ZETA',17X,'V',8X,'VCORR',
C     *      5X,'XINT',7X,'X',9X,'Y')
C Loop over the Beta quadrature point in the nonadiabatic region
         DO 80 IB = IB0,IB1
            BETA = BETAPT(IB)
            DBET = BETA - BETA0
            X = XAB0 + DBET*DXDB
            Y = YAB0 + DBET*DYDB
            R2 = Y/CM2
            R1 = X - CM1*R2
            R3 = R1 + R2
            CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)                     GCL0893 
            PS = V - E
            IF (IC3D .EQ. 2) THEN
               IF (DB.GT.0.0D0) THEN                                    GCL1092
                  VBCORR = VB0 + DBET*DVB/DB
               ELSE
                  VBCORR = VB0 + 0.5D0*DVB                              GCL1092
               END IF
               PS = PS + VBCORR
            END IF
            XINT = 0.0D0                                                GCL1092
            IF (PS .GT. 0.D0) THEN                                      GCL1092
               XINT = DXIDB*SQRT(PS)
               SUM(1) = SUM(1) + WEIGHT(IB,1)*XINT
               T = WEIGHT(IB,2)*XINT
               PERC(1) = PERC(1) + T
               SUM(2) = SUM(2) + T
            END IF
C            WRITE (IDBG,9708) IB,BETA,(1.0D0+BETA)*ZETAMX,V,           GCL1092
C     *        VBCORR,XINT,X,Y
C9708        FORMAT(1X,I5,1P2E10.2,10X,5E10.2)
80       CONTINUE
         IF (IB0.EQ.1) THEN 
            B0 = -1.0D0                                                 GCL1092
         ELSE
            B0 = 0.5D0*(BETAPT(IB0-1)+BETAPT(IB0))                      GCL1092
         END IF
         IF (IB1.EQ.NBETA) THEN 
            B1 = 1.0D0                                                  GCL1092
         ELSE
            B1 = 0.5D0*(BETAPT(IB1+1)+BETAPT(IB1))                      GCL1092
         END IF
         PERC(2) = 50.0D0*(B1 - B0)                                     GCL1092
      ELSE
C Sum up contributions from adiabatic region
         DO 90 IB = 1,NBETA
            SUM(1) = SUM(1) + WEIGHT(IB,1)*XINTG(IB)
            SUM(2) = SUM(2) + WEIGHT(IB,2)*XINTG(IB)
90       CONTINUE
      END IF
      T1 = SQRT(2.0D0*XMU)                                              GCL1092
      THETA(1) = T1*SUM(1)
      THETA(2) = T1*SUM(2)
      IF (SUM(2) .NE. 0.D0) PERC(1) = 100.0D0*PERC(1)/SUM(2)            GCL1092
C      WRITE (IDBG, 9709) SUM,THETA,PERC
C9709  FORMAT(' SUM,THETA,PERC=',1P6E13.5)
6000   FORMAT(/,2X,T5,'Error: In call to LBETAS from LINTG,',
     *          ' ALF,BETA=', F10.5, 1X,F10.5)
      RETURN
      END
