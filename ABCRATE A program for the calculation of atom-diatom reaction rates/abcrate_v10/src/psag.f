!!!*************************************************************
! 文件/File: psag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: psag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PSAG (IC3D, EZ, VM, PE, EP, NDIM)
C
C     PSAG   - compute SAG-type probabilites
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     AITKEN - do Aitken interpolation of V, D, XK, and CUR
C     DERS   - derivatives of morse turning point and zero pt. energy
C              w.r.t. s
C     SAGARG - compute effective mass terms for adiabatic tunneling
C              calcuations
C     SPL1B2 - evaluate spline fit
C     TP     - find turning points in adiabatic barrier
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LRFLC, LRFLC2
      CHARACTER*10 ATYP(4)
      DIMENSION EP(NDIM), PE(2,NDIM,4)
      DIMENSION A(4), STP(10), SUM(2,4), TAB(3), THETA(2,4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /ESPEC/  ESPEC(40), NESPEC
      COMMON /INTER1/ NINT
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PSAGCM/ NPSAG
      COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
      COMMON /QUADTH/ PT2(NQGKDM), WT2(NQGKDM,2), NQ22, NSEG2
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE ATYP, NK1, NKT                                               TCA1097
      DATA ATYP /'MEPSAG', 'MCPSAG', 'CD-SCSAG', 'PASAG'/
      DATA NK1 /4/, NKT /4/
C
      NQMID = (NQ22+1)/2
      NQ2 = (NQ22-1)/2
      NSGMID = (NSEG2+1)/2
      IF (NPSAG .EQ. 0) THEN
         WRITE (6, 600) ATYP(1), NSEG2, NQ2, NSEG2, NQ22
      ELSE
         WRITE (6, 610) NSEG2, NQ2, NSEG2, NQ22
      END IF
      IF (NESPEC .LE. 0) THEN
         NE = NQ12
      ELSE
         NE = NESPEC
         NSEG = 1
      END IF
      IE = 0
C  Loops over energies
      DO 160 ISEG = 1,NSEG
         DO 150 IIE = 1,NE
            IE = IE + 1
            IF (NESPEC .GT. 0) THEN
               E = ESPEC(IE)
               IF (E .GT. VMAX(IC3D)) E = 2.0D0*VMAX(IC3D) - E
            ELSE
               E = EZ + 0.5D0*(2.D0*ISEG - 1.D0 + PT(IIE))*VM/NSEG
            END IF
            NTP =0
C  zero theta
            DO 20 J = 1,NKT
               DO 10 I = 1,2
                  PE(I,IE,J) = 0.5D0
                  IF (J .GT. NK1) GO TO 10
                  THETA(I,J) = 0.D0                                     GCL1092
   10          CONTINUE
   20       CONTINUE
            IF (E .LT. EZ) GO TO 150
C
            IFLG = 0
C  loop over sets of turning points
   30       CONTINUE
               CALL TP(IFLG, E, SL, SR, SN, XT, LRFLC)
               IF (IFLG .EQ. 0) GO TO 120
C
               DO 50 J = 1,NK1
                  DO 40 I = 1,2
                     SUM(I,J) = 0.D0                                    GCL1092
   40             CONTINUE
   50          CONTINUE
               ISINT = 1
               NSG = NSEG2
               IF (LRFLC) NSG = NSGMID
C  loop over quadrature points for S integral
               DO 90 ISEG2 = 1,NSG
                  NQLAST = NQ22
                  LRFLC2 = LRFLC .AND. ISEG2 .EQ. NSGMID .AND.
     *               (-1)**NSEG2 .LT. 0
                  IF( LRFLC2) NQLAST = NQMID
                  DO 80 N = 1,NQLAST
                     S = SL + XT*(2.D0*ISEG2 - 1.D0 + PT2(N))/NSEG2
                     CALL SPL1B2 (NS, SV, AV, BV, CV, DV, S, TAB, 1)
                     PS = TAB(1) - E
                     IF (PS .LE. 0.0D0) GO TO 80
                     CALL AITKEN (ISINT, S, V, D, XK, CP)
                     CALL DERS (ISINT, S, DUDSM, DEDS)
                     CALL SAGARG (PS, S, D, XK, CP, DUDSM, A)
                     T2 = 1.D0                                          GCL1092
                     IF (LRFLC2 .AND. N .EQ. NQMID) T2 = 0.5D0
                     DO 70 J = 1,NK1
                        T = A(J)*T2
                        DO 60 I = 1,2
                           SUM(I,J) = SUM(I,J) + T*WT2(N,I)
   60                   CONTINUE
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
               T = SQRT(2.D0*XMU)*XT/NSEG2
               IF (LSYM) T = 2 .D0*T
C  sum up theta from contribution from different turnign point regions
               DO 110 J = 1,NK1
                  DO 100 I=1,2
                     THETA(I,J) = THETA(I,J) + T*SUM(I,J)
                     IF (THETA(I,J) .GE. 44.D0) THEN
                        PE(I,IE,J) = EXP(-2.0D0*THETA(I,J))
                     ELSE
                        PE(I,IE,J) = 1.D0/(1.D0+EXP(2.D0*THETA(I,J)))
                     END IF
  100             CONTINUE
  110          CONTINUE
               NTP = NTP + 2
               IF (NTP .LE. 10) THEN
                  STP(NTP-1) = SL
                  STP(NTP) = SR
               END IF
               GO TO 30
  120       CONTINUE
            IF (NESPEC .GT. 0) THEN
               IF (ESPEC(IE) .GE. VMAX(IC3D)) THEN
                  E = ESPEC(IE)
                  DO 140 I = 1,2
                     DO 130 J = 1,NK1
                        PE(I,IE,J) = 1.0D0 - PE(I,IE,J)
  130                CONTINUE
  140             CONTINUE
               END IF
            END IF
            E = E*CKCAL
            EP(IE) = E
            NM = MIN(NTP, 10)
            IF (NPSAG .EQ. 0) WRITE (6, 602) E, (PE(I,IE,1), I=1,2),
     *         NTP,(STP(J), J=1,NM)
  150    CONTINUE
  160 CONTINUE
C      WRITE OUT PROBABILITIES
      IF (NPSAG .EQ. 0) THEN
         I0 = 3
            WRITE (6, 604) ATYP(I0)                                     TCA0197
            WRITE (6, 606) NSEG2, NQ2, NSEG2, NQ22                      TCA0197
            DO 180 IE = 1,NE*NSEG
               WRITE (6, 608) EP(IE), (PE(J,IE,I0), J=1,2)              TCA0197
  180       CONTINUE
      END IF
      RETURN
  600 FORMAT (/ T26, A10/ 3X, 'E(kcal)', T20, 'N=', I2, '*', I2, 6X,
     *   'N=', I2, '*', I2, T52, 'Turning Pts')
  602 FORMAT (1X, 1PE13.5, 2X, 2E13.5, 2X, I5, 0P,10F8.3) 
  604 FORMAT (/ 7X, 4(18X, A10))
  606 FORMAT (3X, 'E(kcal)',T20,'N=',I2,'*',I2,6X,'N=',I2,'*',I2)       TCA0197
  608 FORMAT (1X, 1PE13.5, 4(2X, 2E13.5))
  610 FORMAT (' N=', I2, '*', I2, ' and ', I2, '*', I2)
      END
