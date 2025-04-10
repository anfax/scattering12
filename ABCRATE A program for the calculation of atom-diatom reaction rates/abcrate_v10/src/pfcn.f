!!!*************************************************************
! 文件/File: pfcn.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pfcn.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PFCN (IS, IT, AH, LP)
C
C     PFCN   - compute GTS partition functions and free energy
C
C  Called by:
C     EXTRAS - print out extra info about GTS
C     GTST   - compute free energies, CVT and ICVT rates
C
C  Calls:
C     PFCNB  - compute bending part. fcn.
C     PFCNR  - compute rotational part. fcn.
C     PFCNST - compute stretching part. fcn.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LP, LNPR
      CHARACTER*15 AH
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /PFCN1/  QROT, CROT, QB, CB, DELPHI, CSTR, QSTR, QT, G(2)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SIGCOM/ ISIGMA(2)                                         TCA1097
      PARAMETER (NTEMDM=100)                                            TCA0996
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP                                                     TCA0996
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)             TCA0996
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      BET = BETA (IT)
      LNPR = .TRUE.
      IIS = IS
      IF (IS .GT. 0 .AND. IS .LE. NSDM10) THEN
         S = SS(IS)
C  rotational partition function
         CALL PFCNR (BET, XMOMS(IS), ISIGMA(1), QROTCL, QROT)           TCA1097
         CROT = EXP(QROT - QROTCL)
         ALISIG = LOG(DBLE(ISIGMA(1)))                                  TCA1097
C  bending vibrational partition function
         FB = FBS(IS)
         GB = GBS(IS)
         CALL PFCNB (BET, FB, QFBS(IS), GB, QB, CB, DELPHI, LP)
         XLPHI = 1.D35
         IF (FB .GT. 0.D0) THEN
            XLPHI = 0.D0                                                GCL1092
            DELPHI = DELPHI*180.0D0/PI
            IF (QB .LT. 80.D0) THEN
               WB = SQRT(FB*GB)
               XLPHI = TPI*GB/(BET*WB*WB*4.0D0*PI)
            END IF
         END IF
         D = DS(IS)
         XK = XKS(IS)
         V = VS(IS)*CKCAL
      ELSE
         IF (IS .LE. 0) THEN
            IIS = -1
            D = DBC
            XK = XKBC
            V = 0.0D0
            S = -50.0D0
         ELSE
            IIS = -2
            D = DAB
            XK = XKAB
            V = VDELTA*CKCAL
            S = 50.0D0
         END IF
         LNPR = .FALSE.
         FB = 0.0D0
         QROT = -80.D0                                                  GCL1092
         CROT = 0.0D0
         XLPHI = 0.0D0
         DELPHI = 0.0D0
         QB = -80.D0                                                    GCL1092
         CB = 0.0D0
         G(2) = -80.D0                                                  GCL1092
      END IF
C  stretching vibrational partition function
      ISSV = IIS
      CALL PFCNST (BET, IIS, S, D, XK, QSTR, CSTR, LP, LGS(7))
      IF (IIS .EQ. 99999) THEN
         IIS = ISSV
         WRITE (6, 6000) IS, IIS, SS(IIS), DS(IIS), XKS(IIS), VS(IIS)
         STOP 'PFCN 1'
      END IF
      T = CKCAL/BET
      G(1) = V - T*QSTR + CNSTC(IT)
      IF (LNPR) THEN
         QT = QSTR + QROT + ALISIG + QB
         G(2) = V - T*QT + CNST(IT)
      END IF
      IF (LP .AND. PTEMP(IT)) WRITE(6, 600) S, QROT, ISIGMA(1), CROT,   TCA1097
     *   XLPHI, DELPHI, QB, CB, QSTR, CSTR, G, AH                       TCA1097
      IF (LNPR) QROT = QROT + ALISIG
      RETURN
  600 FORMAT (1X, F7.3, 1X, 1PE12.4, 0P,I4, F8.3, 1PE12.4, 0PF9.3,      GCL0992
     *   2(1PE12.4, 0PF7.3), 1P,2E12.4, 2X, A15)                        GCL0992
 6000 FORMAT (2X,T5,'Error: In PFCN there is an error in the call ',
     *            'to PFCNST',/,2X,T12,'IS,IIS,S,D,XK,V=',2I5,1P,4E13.5) 
      END
