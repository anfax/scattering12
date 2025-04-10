!!!*************************************************************
! 文件/File: lag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LAG (IC3D, EZ, VM, PMAX, P0, ESV, ALFSV, ALFSV2, NDIM)
C
C     LAG    - compute LAG probabilities
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     LALFMN - find minimum theta(alf)
C     LALFST - set up zeroth order path for LAG
C     LCSA   - compute cosine factors for LAG calculation
C     PLAG   - compute primitive LAG probabilities from thetas
C     PNORM  - normalizae LAG probabilities
C     PREAD  - read in info needed for LAG calculation
C     PWRITE - write out info needed for LAG calculation
C     SECOND - returns cpu time
C     TITLES - print out a two-line title
C     TP     - find turning points in adiabatic barrier
C     TP2    - find final turning point in adiabatic barrier in product
C              channel
C     VSPLN2 - spline fit of adiabatic potential in product channel
C     WKB    - compute WKB energy levels for stretch
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LRFLC, LVAR, LNONAD
      CHARACTER*2 AA(2)
      DIMENSION PMAX(2,NDIM,8), P0(2,NDIM), ESV(NDIM), ALFSV(2,NDIM),
     * ALFSV2(2,NDIM)
      PARAMETER (ITPMX=20)
      DIMENSION STP(ITPMX)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      DIMENSION THETA(2,NQKPDM,5)
      DIMENSION PSUM(2,NQKPDM,8),PNRM(2,NQKPDM,8)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /ESPEC/  ESPEC(40), NESPEC
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      COMMON /LAG1/ THETMN(2,NQKPDM), THET0(2,NQKPDM), THETH(2,NQKPDM),
     * THET1S(2,NQKPDM), THET1(2,NQKPDM), PE(2,NQKPDM,8),
     * PERCNT(2,2,NQKPDM), CSA0L(2,NQKPDM,4), CSA0R(2,NQKPDM,4),
     * GAMSV(2,NQKPDM), SLG(NQKPDM), SRG(NQKPDM), THETSM(2,4),
     * ALFMN(2), ALFX(10), PERC(2,2), CSAL(2,4), CSAR(2,4),
     * NTPSV(NQKPDM)
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      COMMON /LAGCM3/ NPLAG
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /OPTION/ IOPT(20)
      COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     * NSSP(NSADDM), LRSTRT
      COMMON /SINT/   SINT, PTS(NQGKDM), WTS(NQGKDM), NSS, NINTS
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /THLAG1/ XM1, YM1, X1, Y1, UXM1, UYM1, UX1, UY1, SL, SR,
     * XDIF, XSUM, YDIF, YSUM
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     * DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     * FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     * X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      EQUIVALENCE (THETA, THETMN)
      SAVE AA                                                           TCA1097
      DATA AA /'1D', '3D'/
C
      EVAG = (VM + EZ) * CKCAL
      IEND = 0
      WRITE (6, 600) AA(IC3D)
      CALL TITLES (1, 6, 1)
      NQ3 = (NQ32 - 1) / 2
      NE = NQ12
      IF (NESPEC .GT. 0) NE = NESPEC
      IF (NESPEC .GT. 0) NSEG = 1
      NEMAX = NSEG * NE
C  zero arrays
      DO 20 J = 1,8
         DO 15 IE = 1,NEMAX
            DO 10 I = 1,2
               PSUM(I,IE,J) = 0.D0                                      GCL1092
               PMAX(I,IE,J) = 0.D0                                      GCL1092
10          CONTINUE
15       CONTINUE
20    CONTINUE
      WRITE (6, 604)
      IF (IOPTAU.EQ.0) WRITE (6, 606)
      IF (IOPTAU.NE.0) WRITE (6, 608)
      IF (LGS2(3)) WRITE (6, 610)
C  loop over final-state quantum number
      DO 500 NFINAL = IOPT(2), IOPT(3)
         LNONAD = NFINAL .NE. NSTATE
         IF (NPLAG .EQ. 0) THEN
            IF (IOPT(1) .EQ. 0) THEN
               WRITE (6, 612)
            ELSE
               WRITE (6, 614) NSTATE, NFINAL
            END IF
         ELSE
            WRITE (6, 615) (NSEG3, NQ3, NSEG3, NQ32, I=1,4)
         END IF
C  read in from restart file if restart option set
         IE0 = 0
         IF (LLAGRS .AND. IEND .EQ. 0) THEN
            CALL PREAD(LGS(9), LGS(10), NEMAX, NQKPDM, EZ, VM,          TCA0997
     *       IOPTAU, ESV, ALFSV, ALFSV2, SLG, SRG, THET1, THETA,
     *       NTPSV, PERCNT, CSA0L, CSA0R, NFINAL, IEND, IE0, IC3D)
C IE0 is last energy successfully read in
            IF (IEND .LT. 0) STOP
         END IF
C ALFX and ALFXX store minimum alf from last energy to be used as
C initial guess for the current energy.  ALFX(NTP) is for the NTPth
C segmented path and ALFXX is for the unsegmented path.  They are both
C initialized at 0.1.
         DO 30 I = 1,10
            ALFX(I) = 0.1D0
30       CONTINUE
         ALFXX = 0.1D0
         IF (IOPT(5) .EQ. 0) THEN
            VPNF = VP
         ELSE
            VPNF = VAD(NSMAX,IC3D)
         END IF
C Spline fit adiabatic potential curve in product channel
         IF (LNONAD) CALL VSPLN2 (IC3D, NFINAL, NSMAX, SMAX(IC3D), VPNF)
         IF (NPLAG .EQ.  0) WRITE (6, 616) (NSEG3, NQ3, NSEG3, NQ32,
     *    I=1,4)
C ======================================================================
C loop over energies
         IE = 0
         IEWRT = NEMAX + 1
         DO 400 ISEG = 1,NSEG
            DO 390 IIE = 1,NE
               IE = IE + 1
C don't compute P's for energies already read in
               IF (IE .LE. IE0) THEN
                  ILAGTH = 0
                  DO 50 J = 1,4
                     DO 40 I = 1,2
                        T = EXP(-2.D0 * THETA(I,IE,J))
                        PE(I,J,1) = T/(1.D0+T)
40                   CONTINUE
50                CONTINUE
                  E = ESV(IE)/CKCAL
                  IF (E .LT. VPNF) THEN
                     IF (NPLAG .EQ. 0) WRITE (6, 620) ESV(IE)
                  ELSE IF (E .LT. VAD(1,IC3D) .OR. (.NOT.LNONAD .AND.
     *                E .LT. VAD(NSMAX,IC3D)) .OR. (LNONAD .AND. E
     *                .LT. VV2(NS2))) THEN
                     IF (NPLAG .EQ. 0) WRITE (6, 621) ESV(IE)
                  ELSE
                     IF (NPLAG .EQ. 0)
     *                WRITE (6,618) ESV(IE),((PE(I,J,1),I=1,2),J=2,4),
     *                (PE(I,1,1),I=1,2),ILAGTH
                     IEWRT = MIN(IE,IEWRT)
                  END IF
C    skip to end of loop over energies
                  GO TO 390
               END IF
C     zero arrays
               DO 80 I = 1,2
                  THET1(I,IE) = 0.D0                                    GCL1092
                  DO 60 J = 1,4
                     THETA(I,IE,J) = 0.D0                               GCL1092
                     THETSM(I,J) = 0.D0                                 GCL1092
60                CONTINUE
                  DO 70 J = 1,2
                     PERCNT(J,I,IE) = 0.D0                              GCL1092
70                CONTINUE
80             CONTINUE
               NTP = 0
               ITP = 0
               SLG(IE) = SS(1)
               SRG(IE) = SS(NSMAX)
C
               IERR = 0
C  *********************************************************************
C     compute thetas
C
C     check for special energies
               IF (NESPEC .GT. 0) THEN
                  E = ESPEC(IE)
                  IF (E .GT. VMAX(IC3D)) E = 2.0D0 * VMAX(IC3D) - E
               ELSE
                  E = EZ + 0.5D0 * VM * (2.D0* ISEG - 1.D0 + PT(IIE))
     *             /NSEG
               END IF
C                 WRITE (99, 900) IE, ESV(IE), VAD(1,IC3D)*CKCAL,
C    *               VAD(NSMAX,IC3D)*CKCAL
C 900             FORMAT (' IE,E,VAD(1),VAD(NSMAX)=' / I5, 1P3E15.7)
C
C     if energy is too low skip computation
               IF (E .LT. VPNF) THEN
                  IERR = 1
                  IF (NPLAG .EQ. 0) WRITE (6, 620) E*CKCAL
                  GO TO 320
               END IF
               IF (E .LT. VAD(1,IC3D) .OR. (.NOT.LNONAD .AND. E
     *          .LT. VAD(NSMAX,IC3D)) .OR. (LNONAD .AND. E .LT.
     *          VV2(NS2))) THEN
                  IERR = 1
                  IF (NPLAG .EQ. 0) WRITE (6, 621) E*CKCAL
                  GO TO 320
               END IF
C
               IFLG = 0
               SL = SS(1)
               SR = SS(NSMAX)
               ILAGTH = 0
               LVAR = .FALSE.
C ----------------------------------------------------------------------
C     loop over turning points
90             CONTINUE
               SLO = SL
               SRO = SR
               IF (LVAR) GO TO 150
               CALL TP (IFLG, E, SL, SR, SN, XT, LRFLC)
               IF (LNONAD .AND. (IFLG .EQ. 0 .OR. SR .GT.
     *          SMAX(IC3D))) THEN
                  LVAR = .TRUE.
                  IF (IFLG .EQ. 0) SL = SLO
                  SRSV = SR
                  CALL TP2 (E, SR)
                  IF (SR .LE. SMAX(IC3D)) THEN
                     SR = SRSV
                  ELSE
                     IFLG = 1
                  END IF
               END IF
               IF (IFLG .EQ. 0) GO TO 150
               XT = (SR-SL)*0.5D0
               SN = (SL+SR)*0.5D0
C                    WRITE (99, 901) SL, SR, SN, XT
C 901                FORMAT (' SL,SR,SN,XT=', 1P4E15.7)
               IF (XT .LT. 1.D-16) GO TO 90
               NTP = NTP + 1
               IF (NTP.EQ.1) SLG(IE) = SL
               SRG(IE) = SR
               IF (ITP.LT.ITPMX) THEN
                  STP(ITP+1) = SL
                  STP(ITP+2) = SR
               END IF
               ITP = ITP + 2
C     set up to compute X,Y along alf=0 path
               CALL LALFST (X0, Y0)
C     find alf that gives minimum theta
C                    WRITE (99, 902)
C 902                FORMAT (' ENTER LALFMN')
               ILAGT = 0
               IF (NTP .LE. 10) ALFMN(2) = ALFX(NTP)
               CALL LALFMN (IC3D, E, ALFMN, THETMN(1,IE), THET0(1,
     *          IE), THETH(1,IE), THET1S(1,IE), PERC, ILAGT,
     *          IERR,LNONAD,NFINAL)
               ILAGTH = ILAGTH + ILAGT
C                    WRITE (99, 903) ALFMN, (THETMN(I,IE), I=1,2),
C    *                  PERC, IERR
C 903                FORMAT (' ALFMN =', 1P2E15.7/' THETMN=', 2E15.7/
C    *                  ' PERC=', 4E15.7, ', IERR=', I5)
               IF (IERR .GT. 0) GO TO 320
               IF (NTP .LE. 10) ALFX(NTP) = ALFMN(2)
C     compute cos(angle between MEP and tunneling path) on reactant side
               IF (NTP .EQ. 1) CALL LCSA (CSAL, CSAR, ALFMN, SL,
     *          SR, LRFLC, 1)
               T = 1.0D0
C     for symmetric system increment NTP and set factor to 2
               IF (LSYM .AND. .NOT.LRFLC) THEN
                  T = 2.D0                                              GCL1092
                  NTP = NTP + 1
               END IF
C     save minimum alf
               DO 100 I = 1,2
                  IF(NTP.EQ.1) ALFSV(I,IE) = ALFMN(I)
                  ALFSV2(I,IE) = ALFMN(I)
100            CONTINUE
C     accumulate percent for alf= .5 and 1 paths
               DO 120 J = 1,2
                  DO 110 I = 1,2
                     PERCNT(I,J,IE) = PERCNT(I,J,IE) + T*PERC(I,J)
110               CONTINUE
120            CONTINUE
C     accumulate theta for alf=alfmn and alf=0, 0.5, and 1 paths
               DO 140 J = 1,4
                  DO 130 I = 1,2
                     THETSM(I,J) = THETSM(I,J) + T*THETA(I,IE,J)
130               CONTINUE
140            CONTINUE
               GO TO 90
150            CONTINUE
C                  write (99,9901) e*ckcal, (thetsm(i,4),i=1,2),
C     *               (csal(i,4),csar(i,4),i=1,2)
C9901              format(' e=', 1pe13.5, ' theta(lcg)=', 2e13.5,
C     *               ' csal=', 2e13.5, ' csar=', 2e13.5)
C     end of loop over turning points
C ----------------------------------------------------------------------
               IF (NTP .GT. 0) THEN
C     calculate and write out P(E) for the alf = 0, .5, 1, and minimum
C     paths. If NTP > 1, this is for the segmented path
                  DO 170 J = 1,4
                     DO 160 I = 1,2
                        T = EXP(-2.D0 * THETSM(I,J))
                        PE(I,J,1) = T/(1.D0+T)
160                  CONTINUE
170               CONTINUE
                  IF (NPLAG .EQ. 0) WRITE (6, 618) E*CKCAL,
     *             ((PE(I,J,1), I=1,2), J=2,4), (PE(I,1,1), I=1,2),
     *             ILAGTH,ITP,(STP(I),I=1,MIN(ITP,ITPMX))
C     save cos(angle) on reactant side
                  DO 190 J = 1,4
                     DO 180 I = 1,2
                        CSA0L(I,IE,J) = CSAL(I,J)
180                  CONTINUE
190               CONTINUE
C
                  IF (.NOT.LSYM) THEN
C     compute cos(angle between MEP and tunneling path) on product side
                     CALL LCSA (CSAL, CSAR, ALFMN, SLO, SRO, LRFLC,
     *                2)
C     save cos(angle) on reactant side
                     DO 220 J = 1,4
                        DO 210 I = 1,2
                           CSA0R(I,IE,J) = CSAR(I,J)
210                     CONTINUE
220                  CONTINUE
                  END IF
                  IF (NTP .GT. 1) THEN
C     multiple turning points found
C     percent in potential region needs to be averaged by NTP
                     DO 240 J = 1,2
                        DO 230 I = 1,2
                           PERCNT(I,J,IE) = PERCNT(I,J,IE) / NTP
230                     CONTINUE
240                  CONTINUE
                     SL = SLG(IE)
                     SR = SRG(IE)
                     IF (LSYM) SR = 2.D0*SS(NMID) - SL
                     SN = 0.5D0 * (SL + SR)
                     XT = 0.5D0 * (SR - SL)
C     set up to compute X,Y along unsegmented alf=0 path
                     CALL LALFST (X0, Y0)
C     find alf that minimizes theta for unsegmented path
C                       WRITE (99,902)
                     ILAGTH = 0
                     ALFMN(2) = ALFXX
                     CALL LALFMN (IC3D, E, ALFMN, THETMN(1,IE),
     *                THET0(1,IE), THETH(1,IE), THET1S(1,IE),
     *                PERC, ILAGTH, IERR, LNONAD,NFINAL)
C                       WRITE (99,903) ALFMN, (THETMN(I,IE),I=1,2),
C    *                     PERC, IERR
                     IF (IERR .GT. 0) GO TO 320
                     ALFXX = ALFMN(2)
C     compute and write out P(E) along unsegmented path
                     DO 260 J = 1,4
                        DO 250 I = 1,2
                           T = EXP(-2.D0 * THETA(I,IE,J))
                           PE(I,J,1) = T/(1.D0+T)
250                     CONTINUE
260                  CONTINUE
                     IF (NPLAG .EQ. 0) WRITE (6, 622) ((PE(I,J,1),
     *                I=1,2), J=2,4), (PE(I,1,1), I=1,2), ILAGTH,
     *                2,SL,SR
C     compute cos(angle) for unsegmented path on both sides
                     CALL LCSA(CSAL, CSAR, ALFMN, SL, SR, LSYM, 3)
C                        write (99,9901) e*ckcal, (thet1s(i,ie),i=1,2),
C     *                   (csal(i,4),csar(i,4),i=1,2)
                     THET1(1,IE) = THET1S(1,IE)
                     THET1(2,IE) = THET1S(2,IE)
C     select theta, cos(angles), etc. from segmented and unsegmented
C     paths to give minimum theta
                     DO 310 J = 1,4
                        IF (THETA(2,IE,J) .GE. THETSM(2,J)) THEN
C     minimum theta is for segmented path, reset theta.
                           DO 270 I = 1,2
                              THETA(I,IE,J) = THETSM(I,J)
270                        CONTINUE
                        ELSE
C     minimum theta is for unsegmented path, reset CSA, ALFSV, PERCNT.
                           DO 280 I = 1,2
                              CSA0L(I,IE,J) = CSAL(I,J)
                              CSA0R(I,IE,J) = CSAR(I,J)
280                        CONTINUE
C                             WRITE (99, 9920) J, PERC
C9920                         FORMAT (' J=', I2, ', PERC=', 10X,
C    *                           1P4E13.5)
                           IF (J .EQ. 1) THEN
                              DO 290 I = 1,2
                                 ALFSV2(I,IE) = ALFMN(I)
                                 ALFSV(I,IE) = ALFMN(I)
290                           CONTINUE
                           END IF
                           IF (J .GT. 2) THEN
                              K = J - 2
                              DO 300 I = 1,2
                                 PERCNT(I,K,IE) = PERC(I,K)
300                           CONTINUE
                           END IF
                        END IF
310                  CONTINUE
                  ELSE
                     THET1(1,IE) = THET1S(1,IE)
                     THET1(2,IE) = THET1S(2,IE)
                  END IF
               END IF
320            CONTINUE
C  *********************************************************************
C
C              WRITE (99, 904) IERR
C 904          FORMAT (' AT PWRITE, IERR=', I5)
               IF (IERR .GT. 0) THEN
                  NTP = 0
                  DO 340 J = 1,5
                     DO 330 I = 1,2
                        THETA(I,IE,J) = 100.D0                          GCL1092
330                  CONTINUE
340               CONTINUE
                  DO 350 I = 1,2
                     ALFSV(I,IE) = 0.0D0
                     ALFSV2(I,IE) = 0.0D0
350               CONTINUE
               ELSE
                  IF (NTP .GT. 0) THEN
                     IEWRT = MIN(IEWRT, IE)
                  ELSE
                     IF (NPLAG .EQ.  0) WRITE (6, 623) E*CKCAL
                  END IF
               END IF
               NTPSV(IE) = NTP
               ESV(IE) = E * CKCAL
               IF (LGS2(8))
     *          CALL PWRITE (LGS(9), LGS(10), NEMAX, NQKPDM, EZ, VM,
     *          IOPTAU, ESV, ALFSV, ALFSV2, SLG, SRG, THET1,            TCA0997
     *          THETA, NTPSV, PERCNT, CSA0L, CSA0R, NFINAL, IE, IC3D)
390         CONTINUE
400      CONTINUE
C end of loop over energies
C ======================================================================
         IF (.NOT.LNONAD) THEN
            DO  420 IE = 1,NEMAX
               DO 410 I =  1,2
                  T = EXP(-2.D0*THET0(I,IE))
                  P0(I,IE) = T / (1.D0 + T)
410            CONTINUE
420         CONTINUE
         END IF
         IF (NPLAG .EQ. 0) THEN
            IESKP = IEWRT - 1
            WRITE (6, 624) (NSEG3, NQ3, NSEG3, NQ32, I=1,2)
            IF (IESKP .GT. 0) WRITE (6, 625) ESV(1), ESV(IESKP)
            IF (IEWRT .LE. NEMAX) WRITE (6, 626) (ESV(IE),
     *       ((PERCNT(J,I,IE), J=1,2), I=1,2), (ALFSV (I,IE), I=1,2),
     *       (ALFSV2(I,IE),I=1,2),IE=IEWRT,NEMAX)
         END IF
         CALL PLAG (NSTATE, NFINAL, NEMAX, NQKPDM, ESV, SLG,            TCA0997
     *    SRG, NTPSV, THETA, CSA0L, CSA0R, PE, GAMSV)
         IF (NPLAG .EQ. 0) THEN
            WRITE (6, 628)
            WRITE (6, 638) SINT, NSS, NINTS
            J0 = 2
            WRITE (6, 640) ATYP(NLAG+J0+2), ATYP(NLAG+J0+5), (NSEG3,    TCA0197
     *       NQ3, NSEG3, NQ32, I=1,2)                                   TCA0197
            IF (IESKP .GT. 0) WRITE (6, 625) ESV(1), ESV(IESKP)
            IF (IEWRT .LE. NEMAX) WRITE (6, 642) (ESV(IE),
     *       (PE(I,IE,J0+1), I=1,2), (PE(I,IE,J0+4), I=1,2),            TCA0197
     *       IE=IEWRT,NEMAX)                                            TCA0197
         END IF
         DO 450 J = 1,8
            DO 440 IE = 1,NEMAX
               DO 430 I = 1,2
                  PSUM(I,IE,J) = PSUM(I,IE,J) + PE(I,IE,J)
430            CONTINUE
440         CONTINUE
450      CONTINUE
C  compute normalized probability from sum
         CALL PNORM (NEMAX,NQKPDM,ESV,EVAG,PSUM,PNRM,GAMSV,THET1)       TCA0997
C  for probabilities summed over multiple states, choose max of the
C  current normalized probability and the previous max.
         DO 470 J = 1,8
            DO 460 IE = 1,NEMAX
               DO 455 I = 1,2
                  PMAX(I,IE,J) = MAX(PMAX(I,IE,J),PNRM(I,IE,J))
455            CONTINUE
460         CONTINUE
470      CONTINUE
500   CONTINUE
C  end of loop over final-state quantum number
      IF (NPLAG .EQ. 0) THEN
         WRITE (6, 644)
         WRITE (6, 638) SINT, NSS, NINTS
         J0 = 2
         WRITE (6, 640) ATYP(NLAG+J0+2), ATYP(NLAG+J0+5), (NSEG3, NQ3,  TCA0197
     *    NSEG3, NQ32, I=1,2)                                           TCA0197
         IF (IESKP .GT. 0) WRITE (6, 625) ESV(1), ESV(IESKP)            TCA0197
         IF (IEWRT .LE. NEMAX) WRITE (6, 642) (ESV(IE), (PMAX(I,IE,     TCA0197
     *    J0+1), I=1,2), (PMAX(I,IE,J0+4), I=1,2), IE=IEWRT,NEMAX)      TCA0197
      END IF
      RETURN
600   FORMAT (/,1X, 25('*'), 1X, A2, ' LCG and LAG calculation ',
     *              25('*'))
604   FORMAT (/, 2X, T5, 'The MEP is used for the ALPHA=0 path')
606   FORMAT (/, 2X, T5, 'In the LCG and LAG calculations of GAMMA, ',
     *                   'the vibrational period ',
     *        /, 2X, T5, 'is approximated from the zero-point ',
     *                   'energy level.')
608   FORMAT (/, 1X, T5, 'In the LCG and LAG calculations of GAMMA, ',
     *                   'the vibrational period ',
     *        /, 1X, T5, 'is computed from the derivative of the ',
     *                   'phase integral with respect to E')
610   FORMAT (/, 1X, T5, 'The optimum ALPHA is chosen from the three ',
     *                   'values: 0., 5., and 1.0;',
     *        /, 1X, T5, 'a root search is not performed.')
612   FORMAT (/, 1X, 23('*'), 1X, 'LAG transmission probabilities',
     *           1X, 23('*'))
614   FORMAT (/, 1X, 12('*'), 2X, 'LAG transmission probabilities ',
     *                            'from NSTR = ', I3, ' to', I3, 
     *           1X, 12('*'))
615   FORMAT (/, 1X, T5, 'The number of quadrature points used in ',
     *                   'evaluating P(E), N= ', I2, '*', I2, ' and ', 
     *           I2, '*', I2)
616   FORMAT (/, 1X, T5, 'For energies with more than two turning ',
     *                   'points, there are two lines.',
     *        /, 1X, T5, 'The top line is summed over the ',
     *                   'contributions from the different regions. ',
     *        /, 1X, T5, 'The bottom line is obtained using the two ',
     *                   'most outer turning points.',
     *       //,5X,'E',T12,'ALF=',T22,'0.0',T43,'0.5',T64,'1.0',T83,
     *             'Optimum',T96,'Func',/,1X,'(kcal)',2X,4('N=',I2,'*',
     *           I2,5X,I2,'*',I2,4X),T96,'Calls',3X,'TPS')
618   FORMAT (1X, 1PE9.3, 1X, 4(1X, 2E10.3), 2I4, (T104,0P,4F7.2)) 
620   FORMAT (1X, 1PE12.5, 2X, 'Energy below product asymptote, P(E)',
     *                         ' set to zero')
621   FORMAT (1X, 1PE12.5, 2X, '***** Turning point outside range,',
     *                         ' P(E) set to zero')
622   FORMAT (11X, 4(1X, 1P,2E10.3), I4,T100,I4,0P,2F7.2) 
623   FORMAT (1X, 1PE12.5, 2X, '***** No sets of turning points found,',
     *                         ' P(E) set to zero')
624   FORMAT (/, 1X, T5, 'Percent contribution to THETA from ',
     *                   'nonadiabatic region (%THETA)',
     *        /, 1X, T5, 'and percent of path length in nonadiabatic ',
     *                   'region (%BETA)',
     *        /, 56X,'Optimum ALPHA',/,2X, 'E(kcal)',4X,'ALPHA = 0.5',
     *            6X,'ALPHA = 1.0',9X,'Reactant',7X,'Product',
     *        /,13X,2(3X,'%THETA',3X,'%BETA'),
     *           2X,2(1X,'N=',I2,'*',I2,2X,I2,'*',I2))
625   FORMAT (/, 1X, T5, 'Calculations were not performed for ',
     *                   'energies from', 1PE13.5, ' to', E13.5)
626   FORMAT ((1X,1PE12.5,2(1X,0P,2F8.3),2X,2(1X,2F7.3))) 
628   FORMAT (/, 1X, T5, 'Unnormalized reaction probabilities')
630   FORMAT (/,T28, A10, T56, A10/,3X, 'E(kcal)', 2(5X, 'N=', I2,
     * '*', I2, 6X, 'N=', I2, '*', I2, 3X))
632   FORMAT (1X, 1PE13.5, 2X, 2E13.5, 2X, 2E13.5)
634   FORMAT (/T28, A10, T56, A10, T82, 'GAM**2'/,3X, 'E(kcal)', 
     * 3(5X, 'N=', I2, '*', I2, 6X, 'N=', I2, '*', I2, 3X))
636   FORMAT ((1X, 1PE13.5, 3(2X, 2E13.5)))
638   FORMAT (/, 1X, T5, 'In the following LCG and LAG calculations, ',
     *                   'the integrals over s are ',
     *        /, 1X, T5, 'divided into segements of approximate ',
     *                   'length SINT with',
     *        /, 1X, T5, 'NSS-order Legendre quadrature used in ',
     *                   'each segment.',
     *        /, 1X, T5, 'Interpolation to the grid points is by ', 
     *                   'NINTS-order Aitken interpolation.',
     *       //, 1X, T5, 'SINT = ', F10.5, ',  NSS = ', I5, 
     *                   ',  NINTS = ', I5)
  640  FORMAT (/, T28, A10, T56, A10, /, 3X, 'E(kcal)', 3X, 2(5X, 'N=', TCA0197
     *         I2, '*', I2, 6X, 'N=', I2, '*', I2, 3X))                 TCA0197
642   FORMAT ((1X, 1PE13.5, 2(2X, 2E13.5)))                             TCA0197
644   FORMAT (/, 1X, T5, 'Normalized reaction probabilities')
      END
