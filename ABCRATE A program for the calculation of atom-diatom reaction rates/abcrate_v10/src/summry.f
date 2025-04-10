!!!*************************************************************
! 文件/File: summry.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: summry.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SUMMRY
C
C     SUMMRY - summarize rate constants
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C     EACT   - compute activation energy for temperature range T1 to T2
C     TABL20 - print out table of ratios
C     TABL22 - print out table of rates
C     TITLES - print 2-line title
C     TRANS  - print out transmission coefficients (kappas)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 AH(2)
      CHARACTER*2 AU(2)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /EACT1/   IACT, NT1(10), NT2(10)
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)                            TCA0197
      LOGICAL LLAG, LLAGRS                                              TCA0197
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     * LLAG, LLAGRS                                                     TCA0197
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      PARAMETER (NTEMDM=100)                                            TCA0996
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP                                                     TCA0996
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)             TCA0996
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      PARAMETER (NRATE=5)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      PARAMETER (NDUM=13*NSDM10-NRATE*NKTYP*NTEMDM-4*NKTYP)
      COMMON /BLKRCM/ RATE(NTEMDM,NKTYP,NRATE), EAC(NKTYP), ARRH(NKTYP),GCL1096
     *ENTRP(2,NKTYP)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL1096
      COMMON /GTSTCM/ XK(NTEMDM,NRATE,2), XKPOL(NTEMDM,2),              GCL1096
     *                QSV(NTEMDM,4,5), THCOR(NTEMDM,2), RAT(4)          GCL1096
      CHARACTER*4 AR(NRATE)
      DIMENSION IPX(NKALAB), JPX(NKALAB)                                TCA0197
      SAVE AH, AR, AU                                                   TCA1097
      DATA AH /'1D', '3D'/, AR /'TST', 'CVT', 'ICVT', 'CUS', 'AG'/
      DATA AU /'CM', 'CC'/
C
      NKAP = NKT + 2
C
      IPX(1)=1                                                          TCA0197
      IPX(2)=2                                                          TCA0197
      IPX(3)=3                                                          TCA0197
      IPXC=3                                                            TCA0197
      DO 5 I=3,NKAP                                                     TCA0197
        IF (LPTYP(I)) THEN                                              TCA0197
          IPXC=IPXC+1                                                   TCA0197
          IPX(IPXC)=I                                                   TCA0197
        END IF                                                          TCA0197
    5 CONTINUE                                                          TCA0197
      DO 6 I=1,NKAP                                                     TCA0197
        JPX(I)=IPX(I)                                                   TCA0197
        IF (I .GT. 3) JPX(I)=JPX(I)+1                                   TCA0197
    6 CONTINUE                                                          TCA0197
C
      WRITE (6, 600)
      CALL TITLES (1, 6, 1)
      WRITE (6, 650)                                                    GCL1096
      WRITE (6, 651)                                                    TCA0197
      IF (LLAG) WRITE (6, 652)                                          TCA0197
      WRITE (6, 602)
C  loop over collinear and 3D
      DO 140 IC3D = 1,2
         IF (.NOT.LGS(IC3D+8)) GO TO 140
C        WRITE (6, 6600) ((QSV(IT,J,IC3D),J=1,4),IT=1,NTMAX)
C6600    FORMAT (' SUMMRY: QSV=', 1P4E15.7 : / (13X, 1P4E15.7))
         WRITE (6, 604) AH(IC3D)
C  write out transmission coefficients
         CALL TRANS(IC3D)
         IFLGEQ = 0
         WRITE (6, 606) AU(IC3D)
C  loop over forward and reverse reaction rates
   20    CONTINUE
C  loop over TST, CVT, ICVT, CUS, and AG
            DO 90 JR = 1,NRATE
               ICAG = MIN(JR, 2)
               NX = MIN(9,IPXC)                                         TCA0197
               WRITE (6, 608) (AR(JR), I=1,NX)                          TCA0197
               IF (JR .EQ. 3 .OR. JR .EQ. 5) THEN
                  WRITE (6, 670) (ATYP(IPX(I)),I=4,NX)                  TCA0197
               ELSE
                  WRITE (6, 668) (ATYP(IPX(I)),I=4,NX)                  TCA0197
               END IF
C  loop over temperatures for kappas
               DO 40 IT = 1,NTMAX
                  R1 = XK(IT,JR,IC3D)
                  IF (IFLGEQ .GT. 0) THEN
                     IF (R1 .GT. 0.0D0) THEN
                        R1 = LOG(R1)
                        R1 = EXP(R1-EQUIL(IT,IC3D))
                     ELSE
                        R1 = R1 * EXP(-EQUIL(IT,IC3D))
                     END IF
                  END IF
                  R2 = R1*KAPW(IT,IC3D)
                  R3 = R1
                  IF (JR .NE. 3 .AND. JR. NE. 5)
     *               R3 = R1*KAP(IT,ICAG,IC3D)
                  RATE(IT,1,JR) = R1
                  RATE(IT,2,JR) = R2
                  RATE(IT,3,JR) = R3
                  DO 30 I = 3,NKAP                                      TCA0197
                     RATE(IT,I+1,JR) = R3*KAP(IT,I,IC3D)
   30             CONTINUE
                  WRITE (6, 612) TEMP(IT), (RATE(IT,JPX(I),JR),I=1,NX)  TCA0197
   40          CONTINUE
   90       CONTINUE
C  activation energies
            IF (IACT .GT. 0) THEN
               WRITE (6, 616) AU(IC3D), AU(IC3D)
               DO 120 JR = 1,NRATE
                  ICAG = MIN(JR, 2)
                  NX = MIN(9,IPXC)                                      TCA0197
                  WRITE (6, 645) (AR(JR), I=1,NX)                       TCA0197
                  IF (JR .EQ. 3 .OR. JR .EQ. 5) THEN
                     WRITE (6, 646) (ATYP(IPX(I)),I=4,NX)               TCA0197
                  ELSE
                     WRITE (6, 647) (ATYP(IPX(I)),I=4,NX)               TCA0197
                  END IF
C  first set of kappas
                  CALL EACT(1, NX, NTEMDM, RATE(1,1,JR), EAC, ARRH,     TCA0197
     *               ENTRP, JPX)                                        TCA0197
  120          CONTINUE
            END IF
            IF (IFLGEQ .EQ. 0) THEN
               CALL TABL20(IC3D)
               CALL TABL22(IC3D)
               IF(LSYM) GO TO 130
               WRITE (6, 620) AU(IC3D)
            END IF
            IFLGEQ = IFLGEQ + 1
         IF (IFLGEQ .LT. 2) GO TO 20
  130    CONTINUE
  140 CONTINUE
C
C  Ratios of partition functions
      WRITE (6, 623)
      IF (LGS(9)) THEN
C  collinear
         WRITE (6, 624)
         DO 170 I = 1,2
            I2 = 2*I
            IF (I .EQ. 2) WRITE (6, 626)
            DO 160 IT = 1,NTMAX
               DO 155 J = 1,2
                  RAT(J) = 0.0D0
                  T = QSV(IT,J,1) - QSV(IT,J,I2)
                  IF (T .LT. 80.D0) RAT(J) = EXP(T)
  155          CONTINUE
               P1 = RAT(1)*RAT(2)
               P2 = 1.0D0
               IF (I .EQ. 2) P2 = 1.0D0/THCOR(IT,1)
               P1 = P1*P2
               WRITE (6, 628) TEMP(IT), RAT(1), RAT(2), P2, P1
  160       CONTINUE
  170    CONTINUE
      END IF
      IF (LGS(10)) THEN
C  3D
         WRITE (6, 630)
         DO 200 I = 1,2
            I2 = 2*I+1
            IF (I.EQ.2) WRITE (6, 632)
            DO 190 IT = 1,NTMAX
               DO 180 J = 1,4
                  RAT(J) = 0.0D0
                  T = QSV(IT,J,1) - QSV(IT,J,I2)
                  IF (T .LT. 80.D0) RAT(J) = EXP(T)
  180          CONTINUE
               P1 = RAT(1)*RAT(2)
               P4 = 1.0D0
               IF (I .EQ. 2) P4 = 1.0D0/THCOR(IT,2)
               P2 = P1*RAT(3)*RAT(4)*P4
               WRITE (6, 629) TEMP(IT), RAT(1), RAT(2), P4, P1, RAT(3),
     *            RAT(4), P2
  190       CONTINUE
  200    CONTINUE
      END IF
      RETURN
  600 FORMAT (1H1/ 1X, 57('*'), ' Summary section ', 57('*'))
  602 FORMAT (/ ' Summary of transmission coefficients, rate',
     *   ' constants, Arrhenius parameters, and entropies of',
     *   ' activation')
  604 FORMAT (// 1X, 20('*'), 1X, 'Summary for ', A2, ' reaction')
  606 FORMAT (// 1X, 20('*'), 1X, 'Forward reaction, rates in units of',
     *   1X, A2, '/molec-sec')
  608 FORMAT (/,T10,'T(K)',T19,9(A4,8X))                                TCA0197
  668 FORMAT (1X,T31,'Wigner',T43,'CAG',T55,6(A10,2X))                  TCA0197
  670  FORMAT (1X,T31,'Wigner',T55,6(A10,2X))                           TCA0197
612   FORMAT (6X,F10.3,1P,T17,E12.4,T29,E12.4,T41,E12.4,T53,E12.4,
     *        T65,E12.4,T77,E12.4,T89,E12.4,T101,E12.4,T113,E12.4)
614   FORMAT (/,T10,'T(K)',T19,A4,T31,A4,T43,A4,T55,A4,T67,A4,
     *        T79,A4,T91,A4,T103,A4,T115,A4,
     *        /,T19,A10,T31,A10,T43,A10,T55,A10,T67,A10,T79,A10,
     *        T91,A10,T103,A10,T115,A10)
  616 FORMAT (// 20X, ' Activation energies, Ea (units of kcal);',
     *   ' A factors (units of ', A2, '/molec-sec);', / 20X, ' and',
     *   ' entropy of activation, S0(T) (units of kcal/deg-mol,',
     *   ' standard state is 1 molec/', A2, ')',/)
  620 FORMAT (// 1X, 20('*'), 1X, 'Reverse reaction, rates in units of',
     *   1X, A2, '/molec-sec')
  623 FORMAT (// 1X, 31('*'), ' Ratios of conventional transition',
     *   ' state to CVT partition functions ', 32('*'))
624   FORMAT (/,1X,T5,'Collinear, at CVT GTS (no threshold correction)',
     *       //,1X,T3,'T(K)',T14,'Stretch',T29,'Exp(-BETA*V)',
     *             T44,'Res.Thr.Cor.',T59, 'Product')
626   FORMAT (/,1X,T5,'Collinear, at ICVT GTS',
     *       //,1X,T3,'T(K)',T14,'Stretch',T29,'Exp(-BETA*V)',
     *             T44,'Res.Thr.Cor.',T59, 'Product')
628   FORMAT (1X,F9.2,1P,T11,E15.6,T26,E15.6,T41,E15.6,T56,E15.6,
     *        T71,E15.6,T86,E15.6,T101,E15.6,T116,E15.6)
629   FORMAT (1X,F9.2,1P,T11,7E15.6)                                    TCA0197
630   FORMAT (/,1X,T5,'3D, at CVT GTS (no threshold correction)',
     *        /,T109,'Ratios of K',1H','s',
     *        /,1X,T3,'T(K)',T14,'Stretch',T29,'Exp(-BETA*V)',
     *             T44,'Res.Thr.Cor.',T59, 'Product',T74,'Bend',
     *             T89,'Rot(Sig=1)',T104,'TST/CVT')                     TCA0197
632   FORMAT (/,1X,T5,'3D, at ICVT GTS',/,T109,'Ratios of K',1H','s',
     *        /,1X,T3,'T(K)',T14,'Stretch',T29,'Exp(-BETA*V)',
     *             T44,'Res.Thr.Cor.',T59, 'Product',T74,'Bend',
     *             T89,'Rot(Sig=1)',T104,'TST/CVT')                     TCA0197
642   FORMAT (/,T6,'T(K)',T19,A4,T31,A4,T43,A4,T55,A4,T67,A4,
     *        T79,A4,T91,A4,T103,A4,T115,A4,
     *        /,T19,A10,T31,A10,T43,A10,T55,A10,T67,A10,T79,A10,
     *        T91,A10,T103,A10,T115,A10)
  645 FORMAT (/,T6,'T(K)',T19,9(A4,8X))
  646 FORMAT (1X,T31,'Wigner',T55,6(A10,2X))
  647 FORMAT (1X,T31,'Wigner',T43,'CAG',T55,6(A10,2X))
  650 FORMAT(/,1X,T5,'TST methods:',                                    TCA0197
     *       /,1X,T5,'TST',T10,'-',T12,'transition state theory',       TCA0197
     *       /,T12,'(dividing surface at saddle point)',
     *       /,T5,'CVT',T10,'-',T12,'canonical variational ',
     *            'transition state theory',
     *       /,T5,'ICVT',T10,'-',T12,'improved canonical variational ',
     *            'transition state theory',
     *       /,T5,'CUS',T10,'-',T12,'canonical unified statistical ',
     *            'model ',
     *       /,T5,'AG',T10,'-',T12,'adiabatic ground-state',
     *       /,T12,'(dividing surface at maximum of ground-state ',
     *             'adiabatic potential)')                              TCA1097
  651 FORMAT (/,1X,T5,'Tunneling methods:',                             TCA0197
     *        /,1X,T5,'MEPSAG',T15,'-',T17,'minimum-energy-path ',      TCA0197
     *          'semiclassical adiabatic ground-state method',          TCA1097
     *        /,1X,T5,'CD-SCSAG',T15,'-',T17,'small-curvature ',        TCA1097
     *          'semiclassical adiabatic ground-state method')          TCA1097
  652 FORMAT (  1X,T5,'LAG',T15,'-',T17,'least-action ground-state ',   TCA0197
     *          'method',                                               TCA0197
     *        /,1X,T5,'LCG3',T15,'-',T17,'large-curvature ground-',     TCA0197
     *          'state method',                                         TCA0197
     *        /,1X,T5,'muOMT',T15,'-',T17,'microcanonical optimized ',  TCA0197
     *          'multidimensional tunneling method')                    TCA0997
      END
