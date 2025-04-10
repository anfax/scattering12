!!!*************************************************************
! 文件/File: plag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: plag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PLAG (NREACT, NPROD, NEMAX, NDIM, ESV, SLG, SRG,       TCA0997
     *NTP, THETA, CSAL, CSAR, PE, GAMSV)
C
C     PLAG   - compute primitive LAG probabilities from thetas
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR
C     AITKF2 - checks for "bad" values before calling Aitken
C              interpolator
C     GAUSSQ - compute gaussian quadrature points and weights
C     VIBTAU - compute vibrational period for u motion
C     SPL1B2 - evaluate spline fit
C     VIBTAP - compute vibrational period for u motion for product
C              channel
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LSAME1, LSAME2, LSAME3
      DIMENSION CSAL(2,NDIM,4), CSAR(2,NDIM,4), ESV(NDIM),
     *GAMSV(2,NDIM), PE(2,NDIM,8), SLG(NDIM), SRG(NDIM),
     *THETA(2,NDIM,5), NTP(NDIM)
      DIMENSION ENDPTS(2), F(7), FX(7), GAM(2), GAM0(2), TAB(3),
     *THT(2,4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      DIMENSION SCR2(NQGKDM)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /SINT/   SINT, PTS(NQGKDM), WTS(NQGKDM), NSS, NINTS
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      DO 3 J = 1,6                                                      TCA0997
         DO 2 IE = 1,NEMAX
            DO 1 I = 1,2
               PE(I,IE,J) = 0.D0                                        GCL1092
1           CONTINUE
2        CONTINUE
3     CONTINUE
      ITP = 0
10    CONTINUE
         ITP = ITP + 1
      IF (NTP(ITP) .EQ. 0 .AND. ITP.LE.NEMAX) GO TO 10
      IF (ITP.GT.NEMAX) RETURN
C  compute the RLAG and RLCG P(E)'s
      DO 20 IE = ITP,NEMAX
         DO 15 I = 1,2
            PE(I,IE,1) = EXP(-2.D0 * THETA(I,IE,5))
            PE(I,IE,2) = EXP(-2.D0 * THETA(I,IE,1))
15       CONTINUE
20    CONTINUE
C  numerical integrals over S. LCG3 and LAG methods
      CALL GAUSSQ(1, NSS, 0.D0, 0.D0, 0, ENDPTS, SCR2, PTS, WTS)
C
      DO 90 J = 1,4
         JJ = J + 2                                                     TCA0997
         DO 80 I = 1,2
            PE(I,ITP,JJ) = EXP(-2.D0*THETA(I,ITP,J))
80       CONTINUE
90    CONTINUE
      IF (ITP.EQ.NEMAX) RETURN
      ITPM = ITP + 1
      ISIMIN = ITP - 1
      NP = NINTS + 1
      IF (NP .GT. 7 .OR. NP .GT. NEMAX-ITP+1) THEN
         NP = MIN(NEMAX-ITP+1,7)
         WRITE (6, 6000) NINTS, NP-1
      END IF
      NINTSN = NP - 1
      NN = NP/2
      TMUI = 2.0D0/XMU
C
C  compute contribution to integral over S from SL1 to SL2
      SL2 = SLG(ITPM-1)
      ISINT = 1
C Loop over energy grid points
      DO 200 ITP = ITPM,NEMAX
         SL1 = SL2
         SL2 = SLG(ITP)
         DELS = SL2-SL1
         IF (DELS .LE. 0.D0) GO TO 200
         NSEG = INT(DELS/SINT) + 1                                      GCL0795
         DELS = DELS/NSEG
C Loop over integration points in S
         DO 190 ISEG = 1,NSEG
            DO 180 IS = 1,NSS
               S = SL1 + 0.5D0*DELS*(2.D0*ISEG-1.D0+PTS(IS))
               ISI = MAX(ISIMIN, ITP-NN-1)
               ISI = MIN(NEMAX-NP, ISI)
C Interpolate vib. period and adiabatic energy
               CALL AITKEN (ISINT, S, V, D, XK, CP)
               CALL VIBTAU (0, S, D, XK, PER)
               CALL SPL1B2 (NS, SV, AV, BV, CV, DV, S, TAB, 1)
               TAUI = 1.D0/PER
               T1 = 0.5D0*DELS*TAUI*WTS(IS)
C              WRITE (98, 9820) ISI, S, TAUI, TAB(1)
C9820          FORMAT ( ' ISI,S,TAUI,V=', I5, 1P3E15.7)
C Loop over methods
               DO 140 J = 1,4
                  JJ = J + 2                                            TCA0997
C Interpolate theta from energy grid point to S
                  DO 130 I = 1,2
                     NNP = 0
                     DO 100 K = 1,NP
                        NNP = NNP + 1
                        F(K) = THETA(I,ISI+K,J)
                        IF (F(K) .LE. 0.D0) GO TO 110
  100                CONTINUE
  110                CONTINUE
                     IF (NNP .GT. 1 .AND. ISI+NNP .GE. ITP) THEN
                        NNINT = NNP - 1
                        T2 = AITKF2 (S, F, FX, SLG(ISI+1), NNINT)
                        IF (T2 .LE. 0.D0) THEN
                           T2 = ((S-SL1)*THETA(I,ITP,J) -
     *                        (S-SL2)*THETA(I,ITP-1,J))/(SL2-SL1)
                        END IF
                     ELSE
                        T2 = 0.D0                                       GCL1092
                     END IF
C Interpolate cosine of angle between the reaction path and the
C    tunneling path
                     DO 120 K = 1,NP
                        F(K) = CSAL(I,ISI+K,J)
  120                CONTINUE
                     CSA = AITKF2 (S, F, FX, SLG(ISI+1), NINTSN)
C For the nonadiabatic motion integrand weighted by 1-cos**2, therefore
C    sign of cos doesn't matter
                     CSA = ABS(CSA)
                     CSA = MIN(1.D0, CSA)
                     T = SQRT(1.D0-CSA*CSA)
                     EX = T1*EXP(-T2)
                     THT(I,J) = T*EX
C                    WRITE (98, 9821) I, J, NP, NNP, T2, CSA, EX
C9821                FORMAT (' I,J,NP,NNP,T2,CSA,EX=', 4I5, 1P3E15.7)
  130             CONTINUE
  140          CONTINUE
C loop over energy grids points summing up contribution
               DO 170 IE=ITP,NEMAX
                  T = ESV(IE)/CKCAL - TAB(1)
                  IF (T .GT. 0.D0) THEN
                     VRI = 1.D0/SQRT(TMUI*T)
                     DO 160 J = 1,4
                        JJ = J + 2                                      TCA0997
                        DO 150 I = 1,2
                           PE(I,IE,JJ) = PE(I,IE,JJ) + VRI*THT(I,J)
  150                   CONTINUE
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C     WRITE (98, 9800)
C9800 FORMAT (' CONTRIBUTION FROM REACTANT CHANNEL'/ 4X, 'IE', 5X, 'PE')
C     WRITE (98, 9801) (IE, (PE(2,IE,J+2), J=1,4), IE=2,NEMAX)          TCA0997
C9801 FORMAT (1X, I5, 1P4E13.5)
      IF (LSYM.AND.NREACT.EQ.NPROD) THEN
C        WRITE (98, 9802)
C9802    FORMAT (' LSYM IS TRUE, MULTIPLY BY TWO.')
         DO 220 J = 1,4
            JJ = J + 2                                                  TCA0997
            DO 210 IE = 1,NEMAX
               DO 205 I = 1,2
                  PE(I,IE,JJ) = 2.D0*PE(I,IE,JJ)
  205          CONTINUE
  210       CONTINUE
  220    CONTINUE
      ELSE
C  compute contribution to integral over S from SR1 to SR2
         SR2 = SRG(ITPM-1)
         ISINT = NSMAX
C Loop over energy grid points
         DO 330 ITP = ITPM,NEMAX
            SR1 = SR2
            SR2 = SRG(ITP)
            DELS = SR1-SR2
            IF (DELS .LE. 0.D0) GO TO 330
            NSEG = INT(DELS/SINT) + 1                                   GCL0795
            DELS = DELS/NSEG
C Loop over integration points in S
            DO 320 ISEG = 1,NSEG
               DO 310 IS = 1,NSS
                  S = SR1 - 0.5D0*DELS*(2.D0*ISEG-1.D0+PTS(IS))
                  CALL AITKEN (ISINT, S, V, D, XK, CP)
                  IF (NREACT .NE. NPROD) THEN
C if reactant and product states are not the same compute period and
C    adiabatic energy from product region
                     CALL VIBTAP (NPROD, S, D, XK, PER)
                     NNS2 = NS2 - ISN + 1
                     CALL SPL1B2 (NNS2, SS(ISN), AV2(ISN), BV2(ISN),
     *                  CV2(ISN), DV2(ISN), S, TAB, 1)
                  ELSE
                     CALL VIBTAU (0, S, D, XK, PER)
                     CALL SPL1B2 (NS, SV, AV, BV, CV, DV, S, TAB, 1)
                  END IF
                  TAUI = 1.D0/PER
                  ISI = MAX(ISIMIN, ITP-NN-1)
                  ISI = MIN(NEMAX-NP, ISI)
                  T1 = 0.5D0*TAUI*DELS*WTS(IS)
C                 WRITE (98, 9820) ISI, S, TAUI, TAB(1)
C Loop over methods
                  DO 270 J = 1,4
                     JJ = J + 2                                         TCA0997
                     DO 260 I = 1,2
C Interpolate theta
                        NNP = 0
                        DO 230 K = 1,NP
                           NNP = NNP + 1
                           F(K) = THETA(I,ISI+K,J)
                           IF (F(K) .LE. 0.D0) GO TO 240
  230                   CONTINUE
  240                   CONTINUE
                        IF (NNP .GT. 1 .AND. ISI+NNP .GE. ITP) THEN
                           NNINT = NNP - 1
                           T2 = AITKF2 (S, F, FX, SRG(ISI+1), NNINT)
                           IF (T2 .LE. 0.D0) THEN
                              T2 = ((S-SR1)*THETA(I,ITP,J) -
     *                           (S-SR2)*THETA(I,ITP-1,J))/(SR2-SR1)
                           END IF
                        ELSE
                           T2 = 0.D0                                    GCL1092
                        END IF
C Interpolate cosine of angle between the reaction path and the
C    tunneling path
                        DO 250 K = 1,NP
                           F(K) = CSAR(I,ISI+K,J)
  250                   CONTINUE
                        CSA = AITKF2(S, F, FX, SRG(ISI+1), NINTSN)
C For the nonadiabatic motion integrand weighted by 1-cos**2, therefore
C    sign of cos doesn't matter
                        CSA = ABS(CSA)
                        CSA = MIN(1.D0, CSA)
                        T = SQRT(1.D0-CSA*CSA)
                        EX = T1*EXP(-T2)
                        THT(I,J) = T*EX
C                       WRITE (98, 9821) I, J, NP, NNP, T2, CSA, EX
  260                CONTINUE
  270             CONTINUE
C Loop over energy grid points summing up contributuions
                  DO 300 IE=ITP,NEMAX
                     T = ESV(IE)/CKCAL - TAB(1)
                     IF (T .GT. 0.D0) THEN
                        VRI = 1.D0/SQRT(TMUI*T)
                        DO 290 J = 1,4
                           JJ = J + 2                                   TCA0997
                           DO 280 I = 1,2
                              PE(I,IE,JJ) = PE(I,IE,JJ) + VRI*THT(I,J)
  280                      CONTINUE
  290                   CONTINUE
                     END IF
  300             CONTINUE
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C     WRITE (98, 9803)
C9803 FORMAT (' CONTRIBUTION FROM PRODUCT CHANNEL'/ 4X, 'IE', 5X, 'PE')
C     WRITE (98, 9801) (IE, (PE(2,IE,J+2), J=1,4), IE=2,NEMAX)          TCA0997
      END IF
C     WRITE (98, 9804)
C9804 FORMAT (2X, 'IE', 1X, 'JJ', 2X, 'I', 6X, 'CSA', 9X, 'PV', 11X,
C    *   'P1', 11X, 'PE')
C Add in contribution from adiabatic mechanism
      DO 360 IE = ITPM,NEMAX
         DO 350 J = 1,4
            JJ = J + 2                                                  TCA0997
            DO 340 I = 1,2
               CSA = MIN(1.D0,CSAL(I,IE,J))
               CSA = MAX(0.D0,CSA)
               IF (.NOT. LSYM) THEN
                  T = MIN(1.D0,CSAR(I,IE,J))
                  T = MAX(0.D0,T)
                  CSA = 0.5D0*(T+CSA)
               END IF
               PV = PE(I,IE,JJ)
               PE(I,IE,JJ) = PV*PV
               IF (NREACT .EQ. NPROD) THEN
                  P1 = CSA*EXP(-THETA(I,IE,J))
                  PE(I,IE,JJ) = PE(I,IE,JJ) + P1*P1
               END IF
C              WRITE (98, 9805) IE, JJ, I, CSA, PV, P1, PE(I,IE,JJ)
C9805          FORMAT (1X, 3I3, 1P4E13.5)
  340       CONTINUE
  350    CONTINUE
  360 CONTINUE
      RETURN
  600 FORMAT(2X,T5,'Warning: In PLAG the number of turning points ',
     *             'differ for the ',
     *     /,2X,T14,'first two energies.  GAMMA has been set to one.',/)
  601 FORMAT(2X,T5,'Warning: In PLAG the number of turning points ',
     *             'differ for the ',
     *     /,2X,T14,'last two energies.  GAMMA has been set to zero.',/)
  602 FORMAT(2X,T5,'Warning: In PLAG at E=',1PE13.5,', the number of ',
     *             'turning points differ', 
     *     /,2X,T14,'with the previous and next energies.  GAMMA has ',
     *              'been set to zero.',/)
 6000 FORMAT (2X,T5,'Warning: In PLAG the order of interpolation has ',
     *              'been reset from',I5,' to',I5)
      END
