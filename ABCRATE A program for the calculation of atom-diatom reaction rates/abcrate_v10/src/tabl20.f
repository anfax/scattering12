!!!*************************************************************
! 文件/File: tabl20.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tabl20.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TABL20 (IC3D)
C
C     TABL20 - print out table of ratios
C
C  Called by:
C     SUMMRY - summarize rate constants
C
C  Calls:
C     TITLES - print 2-line title
C     TOUT   - print out table
C     WRITEF - write with FORM= construct
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LMVT, LUS, LWRITE
      CHARACTER*2 AH(2)
      CHARACTER*10 ATYP2(3)
      CHARACTER*18 FORM(3), FORM2
      CHARACTER*102 FORM1
      DIMENSION RRATIO(6)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,
     * LLAG, LLAGRS 
      PARAMETER (NTEMDM=100)                                            TCA0996
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP                                                     TCA0996
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)             TCA0996
      LOGICAL LEQ
      COMMON /RATEQ/ RKQ(NTEMDM,2), RMVT(NTEMDM,2), RUS(NTEMDM,2),
     *LEQ(2)
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      PARAMETER (NRATE=5)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      PARAMETER (NDUM=13*NSDM10-NRATE*NKTYP*NTEMDM-4*NKTYP)
      COMMON /BLKRCM/ RATE(NTEMDM,NKTYP,NRATE), EAC(NKTYP), ARRH(NKTYP),GCL96
     *ENTRP(2,NKTYP)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL1096
C
      SAVE ATYP2, AH, FORM1, FORM, LWRITE                               TCA1097
      DATA ATYP2 /'1', 'WIGNER', 'CAG'/                                 GCL1096
      DATA AH /'1D', '3D'/
      DATA FORM1 /'(10X,  (1H_)/T17,3HT,K,T25,3HTST,T35,3HCVT,T45,3HICV,GCL1096
     *T48,1HT,T55,3HCUS'/
      DATA FORM /'(T36,6HKappa=,A10)','(T43,6HKappa=,A10)',
     *           '(T46,6HKappa=,A10)'/
      DATA LWRITE /.TRUE./
C
      IF (.NOT.LWRITE .OR. .NOT.LEQ(IC3D)) RETURN
      OPEN (UNIT=20, FILE='abc.20', FORM='FORMATTED',                   GCL1096
     *      STATUS='NEW', ERR=100)                                      GCL1096
      WRITE (20, 2000, ERR=5)
      GO TO 6
    5 CONTINUE
         WRITE (6, 6000)
         LWRITE = .FALSE.
         RETURN
    6 CONTINUE
C write out title to table
      CALL TITLES (1, 20, 1)
      WRITE (20, 2001) AH(IC3D)
C check if US of MUVT rates are given
      LUS = .FALSE.
      LMVT = .FALSE.
      DO 10 IT = 1,NTMAX
         LMVT = LMVT .OR. (RMVT(IT,IC3D).NE.0.D0)
         LUS = LUS .OR. (RUS(IT,IC3D).NE.0.D0)
   10 CONTINUE
C get correct format
      IFORM = 1
      IF (LUS .OR. LMVT) IFORM = 2
      IF (LUS .AND. LMVT) IFORM = 3
      IF (IFORM .LE. 1) THEN
         FORM1(6:7) = '55'
         FORM1(71:83) = '/10X,55(1H_))'
      ELSE IF (IFORM .EQ. 2) THEN
         FORM1(6:7) = '65'
         IF (LMVT) FORM1(71:93) = ',T65,3HuVT/10X,65(1H_))'
         IF(LUS) FORM1(71:92) = ',T65,2HUS/10X,65(1H_))'
      ELSE
         FORM1(6:7) = '75'
         FORM1(71:102) = ',T65,3HuVT,T75,2HUS/10X,75(1H_))'
      END IF
      CALL WRITEF (1, 20, FORM1, ATYP)
      FORM2 = FORM(IFORM)
C loop over kappas
      NLOOP = NKT + 3
      DO 40 ITYP = 1,NLOOP
         INDX = ITYP
         IF (ITYP .GT. 3) THEN
            IF (.NOT. LPTYP(ITYP-1)) GOTO 40                            TCA0197
            CALL WRITEF (2, 20, FORM2, ATYP(ITYP-1))
         ELSE
            CALL WRITEF (2, 20, FORM2, ATYP2(ITYP))
         END IF
         DO 30 IT = 1,NTMAX
            R1 = RKQ(IT,IC3D)
            IF (R1 .EQ. 0.D0) GO TO 30
            DO 20 J = 1,4
               RRATIO(J) = RATE(IT,INDX,J)/R1
   20       CONTINUE
            C2 = 1.D0                                                   GCL1092
            IF (ITYP .EQ. 2) C2 = KAPW(IT,IC3D)
            IF (ITYP .GT. 3) C2 = KAP(IT,INDX-1,IC3D)
            II = 4
            IF (LMVT) THEN
               II = II + 1
               RRATIO(II) = RMVT(IT,IC3D)*C2/R1
            END IF
            IF (LUS) THEN
               II = II + 1
               RRATIO(II) = RUS(IT,IC3D)*C2/R1
            END IF
            CALL TOUT(TEMP(IT), RRATIO, II)
   30    CONTINUE
         WRITE (20, 2002)
   40 CONTINUE
      RETURN
C
100   WRITE (6, 900)
      LWRITE = .FALSE.
C
 2000 FORMAT(1H1)
 2001 FORMAT(/// 10X, 'Ratios of Various Thermal Rate Constants to the',
     *   ' Accurate'/ 10X, 'Ones for the ', A2, ' Reaction' ///)
 2002 FORMAT(' ')
900   FORMAT(/,1X,T5,'Error: Cannot open abc.20 in TABL20')
6000  FORMAT(/,2X,T5,'Error: In TABL20 there is a problem writing ',
     *               'the table of rates to ',
     *      /,2X,T12,'the file linked to FORTRAN unit 20',/)
      END
