!!!*************************************************************
! 文件/File: tabl22.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tabl22.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TABL22 (IC3D)
C
C     TABL22 - print out table of rates
C
C  Called by:
C     SUMMRY - summarize rate constants
C
C  Calls:
C     TITLES - print 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LFIRST, LWRITE
      CHARACTER*2 AH(2), AU(2)
      DIMENSION INDX(7)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,
     * LLAG, LLAGRS 
      LOGICAL LGS(10), LGS2(10)                                         GCL1096
      COMMON /LOGIC/ LGS                                                GCL1096
      COMMON /LOGIC2/ LGS2
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
     *   ENTRP(2,NKTYP)                                                 TCA1296
      CHARACTER*4 AR(NRATE)
      SAVE NRP, AR, LFIRST, LWRITE, AH                                  TCA1097
      DATA NRP / 3 /                                                    TCA0197
      DATA AR /'TST', 'CVT', 'ICVT', 'CUS', 'AG'/                       GCL1096
      DATA LFIRST, LWRITE /.TRUE., .TRUE./
      DATA AH /'1D', '3D'/, AU /'cm', 'cc'/
C
      IF (.NOT.LWRITE .OR. .NOT.LGS2(2)) RETURN
      IF (LFIRST) THEN
          OPEN (UNIT = 22, FILE = 'abc.22', FORM = 'FORMATTED',         GCL1096
     *          STATUS = 'NEW', ERR = 100)                              GCL94
         GO TO 6
    5    CONTINUE
            WRITE (6, 6000)
            LWRITE = .FALSE.
            RETURN
    6    CONTINUE
         CALL TITLES (1, 22, 1)
         LFIRST = .FALSE.
         WRITE (22, 2207)                                               GCL1096
         WRITE (22, 2208)                                               TCA0197
         IF (LLAG) WRITE (22, 2209)                                     TCA0197
         IF (LGS(2) .OR. LGS(5)) GO TO 50                               GCL1096
         INDX(1) = 0
         INDX(2) = 3
         INDX(3) = NKT + 1                                              TCA0197
         IF (LLAG) THEN                                                 TCA0197
           NRP = 6                                                      TCA0997
           INDX(4) = NLAG + 4                                           GCL1096
           INDX(5) = NLAG + 7                                           GCL1096
           INDX(6) = NLAG + 8                                           GCL1096
         END IF                                                         TCA0197
      END IF
      IF (LGS(2) .OR. LGS(5)) GO TO 50                                  GCL1096
      WRITE (22, 2201) AH(IC3D), AU(IC3D)
      DO 20 JR = 1,NRATE
         WRITE (22, 2202) (AR(JR),I=1,NRP)                              TCA0197
         WRITE (22,2212) (ATYP(INDX(I)),I=2,NRP)                        TCA0197
         DO 10 IT = 1,NTMAX
            WRITE (22, 2203) TEMP(IT), (RATE(IT,INDX(I)+1,JR),I=1,NRP)  TCA0197
   10    CONTINUE
   20 CONTINUE
      RETURN                                                            GCL1096
C
50    WRITE (22, 2205) AH(IC3D), AU(IC3D)                               GCL1096
      DO 60 IT = 1, NTMAX                                               GCL1096
            WRITE (22, 2206) TEMP(IT), (RATE(IT, 1, I), I = 1, 2)       GCL1096
60    CONTINUE
C
      RETURN
C
100   WRITE (6, 900)
      LWRITE = .FALSE.
C
 2201 FORMAT (//,8X, 'Forward rate constants for ', A2, ' reaction',    TCA1097
     *   ' (units of ', A2, '/molecule-sec)')                           TCA1097
 2202 FORMAT (/,1X,T10,'T(K)',T19,7(A4,8X))                             TCA0197
 2212 FORMAT (1X,T31,6(A10,2X))                                         TCA0197
 2203 FORMAT (6X,F10.3,1P,T17,7E12.4)                                   TCA0197
2205  FORMAT (/,8X, 'Forward rate constants for ', A2, ' reaction',
     *              ' (units of ', A2, '/molecule-sec)',//,
     *        /,1X,T21,'TST',T33,'TST',/,T32,'Wigner')
2206  FORMAT (6X, F10.3, 1P,2E12.4)   
 2207 FORMAT(/,1X,T5,'TST methods:',                                    TCA0197
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
 2208 FORMAT (/,1X,T5,'Tunneling methods:',                             TCA0197
     *        /,1X,T5,'MEPSAG',T15,'-',T17,'minimum-energy-path ',      TCA0197
     *          'semiclassical adiabatic ground-state method',          TCA0197
     *        /,1X,T5,'CD-SCSAG',T15,'-',T17,'small-curvature ',        TCA1097
     *          'semiclassical adiabatic ground-state method')          TCA1097
 2209 FORMAT (  1X,T5,'LAG',T15,'-',T17,'least-action ground-state ',   TCA0197
     *          'method',                                               TCA0197
     *        /,1X,T5,'LCG3',T15,'-',T17,'large-curvature ground-',     TCA0197
     *          'state method',                                         TCA0197
     *        /,1X,T5,'muOMT',T15,'-',T17,'microcanonical optimized ',  TCA0197
     *          'multidimensional tunneling method')                    TCA0997

900   FORMAT(/,1X,T5,'Error: Cannot open abc.22 in TABL22')
6000  FORMAT(/,1X,T5,'Error: In TABL22 there is a problem writing ',
     *               'the table of rate constants ',
     *       /,1X,T12,'to the file linked to FORTRAN unit 22')
      END
