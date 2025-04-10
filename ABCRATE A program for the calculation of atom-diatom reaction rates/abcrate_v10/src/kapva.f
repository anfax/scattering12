!!!*************************************************************
! 文件/File: kapva.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: kapva.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE KAPVA
C
C     KAPVA  - compute kappas
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C
C     BOLTZ  - do numerical integration of P's to get kappas
C     KAPTYP - set up labels and indices for kappas
C     KG1    - set up Kronrod quadrature nodes and weights
C     LAG    - compute LAG probabilities
C     PMUOMT - determine the muOMT probabilities
C     PSAG   - compute SAG-type probabilites
C     PORCPU - returns cpu time
C     TITLES - print out a two-line title
C     VSPLIN - spline fit of adiabatic potential
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION KAPOB(4,6)
      DIMENSION SUM(4), EMAX(2), ALFMX(2), ALFMX2(2), F(3), FX(3)
      DIMENSION IEMAX(2)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /ESPEC/  ESPEC(40), NESPEC
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /OPTION/ IOPT(20)
      COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
      COMMON /QUADTH/ PT2(NQGKDM), WT2(NQGKDM,2), NQ22, NSEG2
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      PARAMETER (NTEMDM=100)
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL1096
      PARAMETER (NPETYP=12)
      DIMENSION PE(2,NQKPDM,NPETYP), P0(2,NQKPDM), ESV(NQKPDM),
     *ALFSV(2,NQKPDM), ALFSV2(2,NQKPDM)
C
      CHARACTER*2 AA(2)
      CHARACTER*10 ATYP2(6)
      DIMENSION INDOB(6), INDLAG(4)
      SAVE AA, ATYP2, INDOB, INDLAG                                     TCA1097
      DATA AA /'1D', '3D'/
      DATA ATYP2 /'MEP', 'MCP', 'SC', 'PA', 2*' '/
      DATA INDOB /1, 2, 3, 4, 0, 0/
      DATA INDLAG /1, 2, 0, 0/
C
      NQ1 = NQ12
      NQ2 = NQ22
      NQ3 = NQ32
      CALL KG1 (NQ1, NQGKDM, PT, WT)
      CALL KG1 (NQ2, NQGKDM, PT2, WT2)
      IF (LLAG) CALL KG1 (NQ3, NQGKDM, PT3, WT3)
      NQ12 = 2*NQ1 + 1
      NQ22 = 2*NQ2 + 1
      NQ32 = 2*NQ3 + 1
C  set up kappa labels and indexes
      CALL KAPTYP
C  loop over 1d and 3d
      DO 200 IC3D = 1,2
         IF (.NOT.LGS(IC3D+8)) GO TO 200
         IF (IOPT(4) .EQ. 0) THEN
            EZ = VR
         ELSE
            EZ = VAD(1,IC3D)
         END IF
         IF (IOPT(5) .EQ. 0) THEN
            EZ = MAX(EZ,VP)
         ELSE
            EZ = MAX(EZ,VAD(NSMAX,IC3D))
         END IF
         WRITE (6, 600) AA(IC3D)
         CALL TITLES (1, 6, 1)
         WRITE (6, 602)
         DO 20 J = 1,NKT
            DO 10 IT = 1,NTMAX
               KAP(IT,J+2,IC3D) = 0.0D0
   10       CONTINUE
   20    CONTINUE
         N = NSADMX(IC3D)
         NSHLF = NSSP(N)
         CALL PORCPU(TIME1)                                             GCL0992
         VM = VMAX(IC3D) - EZ
         IS = ISMAX(IC3D)
         IF (IS .LE. 1 .OR. IS .EQ. NSMAX .OR. VM .LE. 0.0D0) THEN
C  set kappas to 1.0 if no barrier
            DO 50 IT = 1,NTMAX
               DO 30 J = 1,NKT
                  KAP(IT,J+2,IC3D) = 1.0D0
   30          CONTINUE
               DO 40 I =1,NOB
                  KAP(IT,IOB+I-1,IC3D) = 1.0D0
   40          CONTINUE
               CALL PORCPU(TIME1)                                       GCL0992
               KAP(IT,1,IC3D) =
     *            EXP(BETA(IT)*(VAD(NSHLF,IC3D)-VMAX(IC3D)))
   50       CONTINUE
            WRITE (6, 604)
         ELSE
C  spline fit adiabatic potential curve
            CALL VSPLIN (NSMAX, SS, SMAX(IC3D), VMAX(IC3D), VAD(1,IC3D))
C  compute adiabatic transmission probabilities
            CALL PSAG (IC3D, EZ, VM, PE, ESV, NQKPDM)
            CALL PORCPU (TIME2)                                         GCL0992
            TIME = TIME2-TIME1
            WRITE (6, 622) TIME
            IF (LLAG) THEN
C  compute LAG probabilities
               TIME1 = TIME2
               CALL LAG (IC3D, EZ, VM, PE(1,1,NLAG), P0, ESV,
     *            ALFSV, ALFSV2, NQKPDM)
               CALL PORCPU (TIME2)                                      GCL0992
               TIME = TIME2-TIME1
               WRITE (6, 622) TIME
C   compute the muOMT probabilities
C   The muOMT probabilities will be stored in PE(2,NQKPDM,11)
            CALL PMUOMT(ESV, PE, NQKPDM, NPETYP)                        TCA0997
            END IF
            IF (NESPEC .LE. 0) THEN
C  Boltzmann average
               WRITE (6, 606) AA(IC3D)
               CALL TITLES (1, 6, 1)
               WRITE (6, 608)                                           GCL0494
               WRITE (6, 609) (NSEG2, NQ2, NSEG2, NQ22, I=1,3), NSEG,
     *            NQ1, NSEG, NQ1, NSEG, NQ12, NSEG, NQ12
               IF (LLAG) WRITE (6, 610) (NSEG3, NQ3, NSEG3, NQ32,I=1,5)
               DO 170 IT = 1,NTMAX
                  BET = BETA(IT)
                  WRITE (6, 612) TEMP(IT)
                  DO 60 I = 1,NKT
                     KAP(IT,I,IC3D) = 1.0D0
   60             CONTINUE
C  overbarrier contributions
                  DO 70 I = 1,NOB
                     II = INDOB(I)
                     CALL BOLTZ (BET, VM, PE(1,1,II), KAPOB(1,I), 1,
     *                  ESV, EMAX, IEMAX, VMAX(IC3D))
                     KAP(IT,IOB+I-1,IC3D) = KAPOB(4,I)
C    Turn off the printing of the MCP and PA information
                     IF (I .NE. 2 .AND. I .NE. 4)                       TCA0197
     *               WRITE (6, 614) ATYP2(I), (KAPOB(J,I),J=1,4), EMAX
   70             CONTINUE
C
C  tunneling contributions, loop over the types of probs (NKP)
C  The variable NKPINT will be used to indicate the number of tunneling
C  probabilities for which integration is necessary.  If no LCT then 
C  NKPINT = NKP, but if LCT then NKPINT = NKP - 1 because the COMT 
C  method does not need to be integrated.
                  WRITE (6, 616)
                  IF (LLAG) THEN                                        GCL93
                      NKPINT = NKP - 1                                  GCL93
                  ELSE                                                  GCL93
                      NKPINT = NKP                                      GCL93
                  ENDIF                                                 GCL93
C
                  DO 160 J = 1, NKPINT                                  GCL93
                     IF (J .GE. NLAG) THEN
                        IF (.NOT.LLAG) GO TO 160
                        IF (J .LE. NLAG) THEN
                           DO 80 I = 1,4
                              KAP(IT,NLAG+I+1,IC3D) = 1.0D0
   80                      CONTINUE
C  LAG alpha=0 method
                           CALL BOLTZ (BET, VM, P0, SUM, 0, ESV,
     *                        EMAX, IEMAX, VMAX(IC3D))
                           DO 90 I = 1,4
                              SUM(I) = SUM(I) + KAPOB(I,1)
   90                      CONTINUE
                           WRITE (6, 618) SUM, EMAX
                        END IF
                     END IF
C
                     CALL BOLTZ (BET, VM, PE(1,1,J), SUM, 0, ESV,
     *                  EMAX, IEMAX, VMAX(IC3D))
C  add in MEP overbarrier contribution
                     DO 100 I = 1,4
                        SUM(I) = SUM(I) + KAPOB(I,1)
  100                CONTINUE
                     KAP(IT,J+2,IC3D) = SUM(4)
                     IF (J .GE. NLAG .AND. J .EQ. NLAG+INDLAG(1) .OR.
     *                  J .EQ. NLAG+INDLAG(2)) THEN
C  Compute alphas at Emax for LAG and RLAG
                        DO 130 I = 1,2
                           ALFMX(I) = 0.D0                              GCL1092
                           ALFMX2(I) = 0.D0                             GCL1092
                           IF (EMAX(I) .LE. 0.D0) GO TO 130
                           II = IEMAX(I) - 1
                           DO 110 K = 1,3
                              F(K) = ALFSV(I,II+K-1)
  110                      CONTINUE
                           ALFMX(I) = AITKF2(EMAX(I), F, FX, ESV(II), 2)
                           DO 120 K = 1,3
                              F(K) = ALFSV2(I,II+K-1)
  120                      CONTINUE
                           ALFMX2(I) = AITKF2(EMAX(I), F, FX, ESV(II),
     *                        2)
  130                   CONTINUE
C   Turn off MCP and PA information
                        IF (LPTYP(J+2))                                 TCA0197
     *                  WRITE (6, 614) ATYP(J+2), SUM, EMAX, ALFMX,
     *                     ALFMX2
                     ELSE
C   Turn off MCP and PA information
                        IF (LPTYP(J+2))                                 TCA0197
     *                  WRITE (6, 614) ATYP(J+2), SUM, EMAX
                     END IF
                     DO 150 K = 2,NOB
C  test to see if other types of overbarrier should be used
                        II = INDOB(K)
                        IF (J .NE. II) GO TO 150
                        DO 140 I = 1,4
                           SUM(I) = SUM(I) - KAPOB(I,1) + KAPOB(I,K)
  140                   CONTINUE
                        KAP(IT,NKP+K+1,IC3D) = SUM(4)
                        IF (LPTYP(NKP+K+1))                             TCA0197
     *                    WRITE (6, 614) ATYP(NKP+K+1), SUM, EMAX
  150                CONTINUE
  160             CONTINUE
C  compute the COMT transmission coefficient
C  The kappa for the LCG3 method is in the NLAG+7 position and the kappa for
C  the SCSAG method is in the NKT+1 position in the array KAP.
C  NLAG and NKT are defined in the subprogarm KAPTYP.  
C
                  IF (LLAG) THEN                                        TCA0997
                      KAP(IT,NLAG+9,IC3D) = MAX(KAP(IT,NLAG+7,IC3D),
     1                                          KAP(IT,NKT+1,IC3D))
                  ENDIF
C  compute TST/CAG kappa
                  VDIF = VAD(NSHLF,IC3D) - VMAX(IC3D)
                  KAP(IT,1,IC3D) = EXP(BET*VDIF)
                  WRITE (6, 620) ATYP(1), KAP(IT,1,IC3D)                GCL93
  170          CONTINUE
            END IF
         END IF
  200 CONTINUE
      RETURN
  600 FORMAT (/,1X, 12('*'), 1X, A2, ' semiclassical adiabatic',
     *              ' transmission probabilities ', 12(1H*))
  602 FORMAT (/ 1X, 'N = number of quadrature points used in',
     *    ' evaluating P(E)')
  604 FORMAT (/,2X,T5,'Warning: No barrier therefore KAPS set to one.')
  606 FORMAT (/,1X, 24('*'),1X,A2, ' transmission coefficients ',
     *              24('*'))
  608 FORMAT (/, 1X, T5, 'KAP', T14, '-', T16, 
     *                   'transmission coefficient',
     *        /, 1X, T5, 'EMAX', T14, '-', T16,
     *                   'location of the maximum of P(E)EXP(-BETA*E)',
     *        /, 1X, T5, 'ALPHAR', T14, '-', T16, 
     *                   'reactant value of the optimum ALPHA at EMAX',
     *        /, 1X, T5, 'ALPHAP', T14, '-', T16, 
     *                   'product value of the optimum ALPHA at EMAX',
     *        /, 1X, T8, 'TST/CAG',
     *        /, 1X, T5, 'KAP', T14, '-', T16, 
     *                   'the classical transmission coefficient for',
     *        /, 1X, T16, 'correcting the classical threshold of ',
     *                    'conventional TST.')
609   FORMAT (/, 1X, T5, 'All tunneling corrections given below ',
     *                   'are normalized',
     *        /, 1X, T5, 'to the maximum in the ground state ',
     *                   'adiabatic curve.',
     *        /, 1X, T5, 'To obtain the proper KAP for conventional ',
     *                   'or CVT TST,',
     *        /, 1X, T5, 'multiply by the corresponding classical ', 
     *                   'transmission coefficient',
     *       //, 1X, T5, 'N = the number of quadrature points used ',
     *                   'in evaluating P(E)',
     *        /, 1X, T5, 'M = the number of quadrature points used ',
     *                   'in the Boltzmann average',
     *       // 35X, 'KAP', T71, 'EMAX(kcal)', T96, 'ALPHAR', T118, 
     *      'ALPHAP'/8X, 4(6X, 2HN=, I2, 1H*, I2), 4X, 2(4X, 2HN=, I2, 
     *      1H*, I2)/ 8X, 4(6X, 2HM=, I2, 1H*, I2))
  610 FORMAT (1X, 'For LAG', 4(6X, 2HN=, I2, 1H*, I2), 4X, 2(4X, 2HN=,
     *   I2, 1H*, I2), 1X, 2(1X, 2(3X, 2HN=, I2, 1H*, I2)))
  612 FORMAT (/ 1X, 'Temp=', F9.2/ 1X, 'Over barrier')
  614 FORMAT (1X, A10, 1P,4E13.5, 2X, 2E11.3, 1X, 2E10.2, 1X, 2E10.2) 
  616 FORMAT (1X, 'Tunneling + over barrier')
  618 FORMAT (1X, 'ALPHA=0', 3X, 1P,4E13.5, 2X, 2E11.3)              
  620 FORMAT (1X,A10, T51, 1PE13.5)                                 
  622 FORMAT (/,1X,T5,'Time for the calculation of the ',
     *                'transmission probabilities', F12.2, ' sec')
      END
