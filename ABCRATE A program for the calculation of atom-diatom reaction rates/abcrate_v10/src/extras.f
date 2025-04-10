!!!*************************************************************
! 文件/File: extras.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: extras.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE EXTRAS (ISMAX, IT)
C
C     EXTRAS - print out extra info about GTS
C
C  Called by:
C     GTST   - compute free energies, CVT and ICVT rates
C
C  Calls:
C     PFCN   - compute GTS partition functions and free energy
C     THRCOR - compute threshold corrections for ICVT
C     VPHI   - compute bending potential as function of angle
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*15 AH(5), AB
      DIMENSION QTOT(2), ISMAX(4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /PFCN1/  QROT, CROT, QB, CB, DELPHI, CSTR, QSTR, QT, G(2)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
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
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL1096
      COMMON /GTSTCM/ XK(NTEMDM,NRATE,2), XKPOL(NTEMDM,2),              GCL1096
     *                QSV(NTEMDM,4,5), THCOR(NTEMDM,2), RAT(4)          GCL1096
      DIMENSION RAT2(3,5), VB(3,5), SCR(5,5), SCR2(10,5)                TCA1097
      EQUIVALENCE (QTOT, QSTR)
      SAVE AH                                                           TCA1097
      DATA AH/'SADDLE POINT', 'CVT(COLLINEAR)', 'CVT(3D)',
     *   'ICVT(COLLINEAR)', 'ICVT(3D)'/
C
C     WRITE (6, 6601) IT, ISMAX
C6601 FORMAT (' EXTRAS: IT, ISMAX=', 10I5)
      BET = BETA(IT)
      RT = CKCAL/BET
C  zero storeage arrays
      DO 30 I = 1,5
         DO 10  J = 1,3
            QSV(IT,J,I) = 0.D0                                          GCL1092
            RAT2(J,I) = 0.D0                                            GCL1092
            VB(J,I) = 0.D0                                              GCL1092
            SCR(J,I) = 0.D0                                             GCL1092
            SCR2(J,I) = 0.D0                                            GCL1092
   10    CONTINUE
         QSV(IT,4,I) = 0.D0                                             GCL1092
         SCR(4,I) = 0.D0                                                GCL1092
         SCR(5,I) = 0.D0                                                GCL1092
         DO 20 J = 6,11
            SCR2(J,I) = 0.D0                                            GCL1092
   20    CONTINUE
   30 CONTINUE
C  Loop over types of transition states.  I=1 is saddle point; I=2,4 are
C     collinear CVT,ICVT; and I=3,5 are 3d CVT,ICVT.
      DO 40 I = 1,5
         IF (I .GT. 1) THEN
            II = MOD (I, 2) + 9
            IF (.NOT.LGS(II)) GO TO 40
         END IF
         IM2 = MOD(I,2) + 1
         IF (I .GT. 1) THEN
            IS = ISMAX(I-1)
         ELSE
            IS = 0
            IF (NSAD .NE. 0) THEN
               N = NSADMX(2)
               IS = NSSP(N)
            END IF
         END IF
C  Evaluate partition functions
         CALL PFCN(IS, IT, AB, .FALSE.)
         IF (I .GE. 4) THEN
            CALL THRCOR(IM2, IS, BET, TCOR, QTOT(IM2))
            IF (TCOR .LE. 0.D0) TCOR = 1.0D0
            THCOR(IT, IM2) = TCOR
            G(IM2) = G(IM2) - RT*LOG(TCOR)
         END IF
C        WRITE (6, 6602) I, IS, QROT, QB, QSTR
C6602    FORMAT (' EXTRAS:, I,IS,QROT,QB,QSTR=', 2I5, 1P3E15.7)
         IF (IS .LE. 0) THEN
            S = -50.D0                                                  GCL1092
            R1 = 50.D0                                                  GCL1092
            R2 = R2ASY
            Y = CM2*R2
            X = R1+CM1*R2
            V = 0.D0                                                    GCL1092
            VT = 0.D0                                                   GCL1092
            WSTR = 4.D0*DBC*CCM/XKBC
         ELSE IF (IS .GT. NSDM10) THEN
            S = 50.D0                                                   GCL1092
            R1 = R1ASY
            R2 = 50.D0                                                  GCL1092
            Y = CM2*R2
            X = R1+CM1*R2
            V = VDELTA
            VT = VDELTA*CKCAL
            WSTR = 4.D0*DAB*CCM/XKAB
         ELSE
            V = VS(IS)
            VT = V*CKCAL
            R2 = Y2(IS) / CM2
            R1 = X2(IS) - CM1*R2
            PHI = PI*(1.D0-DELPHI/360.D0)
            CALL VPHI(R1, R2, PHI, V)
            VB(1,I) = V*CKCAL
            T = PHI - PI
            T = T * T
            FB = CKCAL*FBS(IS)
            QFB = QFBS(IS)*CKCAL
            VB(2,I) = 0.5D0*FB*T
            VB(3,I) = VB(2,I) + QFB*T*T/24.D0                           GCL1092
            SCR(1,I) = X2(IS)
            SCR(2,I) = R2+ CM1P*R1
            SCR(3,I) = FB
            SCR(4,I) = QFB
            SCR(5,I) = PHI*180.D0/PI
            FB = FBS(IS)
            GB = GBS(IS)
            WB = 0.D0                                                   GCL1092
            IF (FB .GT. 0.D0) WB = SQRT(FB*GB)*CCM
            WSTR = 4.D0*DS(IS)*CCM/XKS(IS)
            S = SS(IS)
            X = X2(IS)
            Y = Y2(IS)
            V = VS(IS)
C           WRITE (6, 6603) IT, I, QROT, QB
C6603       FORMAT (' IT,I,QROT,QB=', 2I5, 1P2E15.7)
            QSV(IT,4,I) = QROT
            QSV(IT,3,I) = QB
            RAT2(1,I) = 0.D0                                            GCL1092
            IF (CB .NE. 0.D0) RAT2(1,I) = 1.D0/CB
            RAT2(3,I) = 1.D0/CROT                                       GCL1096
            SCR2(8,I) = WB
            SCR2(9,I) = GB
            SCR2(10,I) = XMOMS(IS)
         END IF
         QSV(IT,2,I) = -BET*V
         QSV(IT,1,I) = QSTR
         RAT2(2,I) = 1.D0/CSTR                                          GCL1096
         SCR2(1,I) = S
         SCR2(2,I) = R1
         SCR2(3,I) = R2
         SCR2(4,I) = X
         SCR2(5,I) = Y
         SCR2(6,I) = VT
         SCR2(7,I) = WSTR
C        WRITE (6, 6600) I, (QSV(IT,J,I),J=1,4)
C6600    FORMAT (' EXTRAS: I, QSV=', I5, 1P4E15.7)
   40 CONTINUE
      IF (PTEMP(IT)) WRITE (6, 614)                                     TCA0996
      DO 50 I = 1,5
         IF (I .GT. 1) THEN
            II = MOD (I, 2) + 9
            IF (.NOT.LGS(II)) GO TO 50
         END IF
         IF (PTEMP(IT)) WRITE (6, 615) (SCR2(J,I),J=1,11), AH(I)        TCA0996
   50 CONTINUE
      IF (PTEMP(IT)) WRITE (6, 616)                                     TCA0996
      DO 60 I = 1,5
         IF (I .GT. 1) THEN
            II = MOD (I, 2) + 9
            IF (.NOT.LGS(II)) GO TO 60
         END IF
         IF (PTEMP(IT)) WRITE (6, 617) AH(I), (SCR(J,I),J=1,5),         TCA0996
     *      (VB(J,I),J=1,3)
   60 CONTINUE
      IF (PTEMP(IT)) WRITE (6, 618)                                     TCA0996
      DO 70 I = 1,5
         IF (I .GT. 1) THEN
            II = MOD (I, 2) + 9
            IF (.NOT.LGS(II)) GO TO 70
         END IF
         IF (PTEMP(IT)) WRITE (6, 619) AH(I), (RAT2(J,I),J=1,3),        TCA0996
     *      SCR2(1,I)
   70 CONTINUE
C
C  Print out table of transition state properties
      DO 80 IC3D = 1,2
         IF (LGS(IC3D+8)) CALL TABL21(TEMP(IT), ISMAX(IC3D), 2, IC3D)
   80 CONTINUE
      RETURN
  614 FORMAT (/ T101, 'Moment', T113, 'Rot'/ T6, 's', T14, 'RAB', T22,
     *   'RAB', T32, 'X', T40, 'Y', T50, 'V', T60, 'Hbar W(str)', T73,
     *   'Hbar W(bnd)', T88, 'GPHI', T99, 'of inertia', T113, 'Sig'/
     *   T13, '(bohr)', T21, '(bohr)', T30, '(bohr)', T38, '(bohr)',
     *   T46, '(kcal)', T61, '(cm**-1)', T74, '(cm**-1)', T87,
     *   '(a.u.)', T101, '(a.u.)')
  615 FORMAT (1X, F8.4, 1X, 2F8.4, 1X, 2F8.4, 1X, 1P,5E13.5, 0PF6.1,    GCL0992
     *   2X, A15)                                                       GCL0992
  616 FORMAT (/ T97, 'Bending potentials'/ T24, 'RA,BC', T36, 'RAB,C',
     *   T52, 'FPHI', T65, 'APHI', T77, 'PHI', T88, 'V(PHI)-V(PI)',
     *   T104, 'VHARM(PHI)', T118, 'VQUARTIC(PHI)'/ T24, '(bohr)', T36, TCA0197
     *   '(bohr)', T51, '(kcal)', T64,'(kcal)',T76, '(deg)', T91, 
     *   '(kcal)',T104,'(kcal)',T112,'(kcal)')
  617 FORMAT (5X, A15, 0P,2F12.6, 2X, 1P,2E13.5, 0PF12.6, 2X, 1P,3E14.5)GCL0992
  618 FORMAT (/,T43, 'Ratios of partition functions at GTST',           GCL1096
     *        /,T32, 'Bend', T51, 'Stretch', T70, 'Rotation', T90,      GCL1096
     *        's'/ T27, '(Harm/Quart)', T49, '(Harm/Morse)',            GCL1096
     *        T67, '(Class/Quant)')
  619 FORMAT (5X, A15, 1P,3E20.10, 5X, 0PF8.4)                          GCL0992
      END
