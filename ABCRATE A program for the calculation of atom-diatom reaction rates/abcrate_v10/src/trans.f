!!!*************************************************************
! 文件/File: trans.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: trans.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE TRANS (IC3D)
C
C     TRANS  - print out transmission coefficients (kappas)
C
C  Called by:
C     SUMMRY - summarize rate constants
C
C  Calls:
C     FRMT1  - computes format for printing kappas
C     FRMT2  - computes format for printing ratios of rates
C     WRITEG - write with FORM= construct
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LERR
      DOUBLE PRECISION KOB
      CHARACTER*4 AR(4)
      CHARACTER*9 AH(2)
      CHARACTER*200 FORM
      DIMENSION OUTPUT(20)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      COMMON /STATE/  TNP1, LSTATE, NSTATE
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
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      DIMENSION IPX(NKATYP)                                             TCA0197
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL96
      SAVE AH, AR, NRATE                                                TCA1097
      DATA AH /'COLLINEAR', '3D'/, AR /'TST', 'CVT', 'ICVT', 'CUS'/
      DATA NRATE /4/
C
      NKAP = NKT+2
      NKZ = NKP+2
      KDA = IOB+5
C
      IPXC=0                                                            TCA0197
      DO 5 I=1,NKAP                                                     TCA0197
        IF (LPTYP(I)) THEN                                              TCA0197
          IPXC=IPXC+1                                                   TCA0197
          IPX(IPXC)=I                                                   TCA0197
        END IF                                                          TCA0197
    5 CONTINUE                                                          TCA0197
C
      N = MIN(6, IPXC)                                                  TCA0197
      WRITE (6, 602) (ATYP(IPX(I)),I=1,N)                               TCA0197
      WRITE (6, 603)
      DO 30 IT = 1,NTMAX
         FORM = '(1X,F9.3'
         IFORM = 12
         OUTPUT(1) = TEMP(IT)
         OUTPUT(2) = KAP(IT,1,IC3D)
         OUTPUT(3) = KAP(IT,2,IC3D)
         OUTPUT(4) = KAPW(IT,IC3D)
         KOB = KAP(IT, IOB, IC3D)
         OUTPUT(5) = KOB
         DO 10 I = 2,5
            CALL FRMT1 (FORM, IFORM, OUTPUT(I))
   10    CONTINUE
         II = 5
         DO 20 I = 3,N
            II = II + 1
            OUTPUT(II) = KAP(IT,IPX(I),IC3D)                            TCA0197
            CALL FRMT1 (FORM, IFORM, OUTPUT(II))
            II = II + 1
            OUTPUT(II) = 0.D0                                           GCL1092
            IF (OUTPUT(II-1) .GT. 0.0D0) OUTPUT(II) = (1.0D0 -
     *         KOB/OUTPUT(II-1))*100.0D0
            CALL FRMT2(FORM, IFORM, OUTPUT(II))
   20    CONTINUE
         FORM(IFORM:IFORM) = ')'
         CALL WRITEG (FORM, LERR, OUTPUT, II)
         IF (LERR) THEN
            WRITE (6, 6000)
            RETURN
         END IF
   30 CONTINUE
   35 CONTINUE
         I0 = N + 1
         WRITE (6, 600)
         N = MIN(IPXC, I0+3)                                            TCA0197
         IF (I0 .GT. N) GO TO 60
         WRITE (6, 604) (ATYP(IPX(I)),I=I0,N)                           TCA0197
         WRITE (6, 603)
         DO 50 IT = 1,NTMAX
            KOB = KAP(IT,IOB,IC3D)
            FORM = '(1X,F9.3,44X'
            IFORM = 13
            OUTPUT(1) = TEMP(IT)
            II = 1
            DO 40 I = I0,N
               IF (I .GT. NKZ) THEN
                  K = IOB + IPX(I) - NKZ                                TCA0197
                  IF (K .GT. KDA) K = KDA
                  KOB = KAP(IT,K,IC3D)
               END IF
               II = II + 1
               OUTPUT(II) = KAP(IT,IPX(I),IC3D)
               CALL FRMT1 (FORM, IFORM, OUTPUT(II))
               II = II + 1
               OUTPUT(II) = 0.D0                                        GCL1092
               IF (OUTPUT(II-1) .GT. 0.0D0) OUTPUT(II) =
     *            (1.D0-KOB/OUTPUT(II-1))*100.0D0
               CALL FRMT2 (FORM, IFORM, OUTPUT(II))
   40       CONTINUE
            FORM(IFORM:IFORM) = ')'
            CALL WRITEG (FORM, LERR, OUTPUT, II)
            IF (LERR) THEN
               WRITE (6, 6000)
               RETURN
            END IF
   50    CONTINUE
      GO TO 35
   60 CONTINUE
      RETURN
  600 FORMAT (' ')
  602 FORMAT (// 20X, 'Transmission coefficients'// T6, 'T,K', T16,
     *   2A10, 2X, 'Wigner', 2X, 'MEP over bar.', 4(4X, A10, 4X))
  603 FORMAT (56X, 4(4X, 'Kappa', 2X, '%Tun', 3X) )
  604 FORMAT (54X, 4(8X, A10))
 6000 FORMAT (2X,T5,'Error with format in subroutine TRANS',/,1X, 
     *        A120,/,1X,A80)
      END
