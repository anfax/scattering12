!!!*************************************************************
! 文件/File: thrcor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: thrcor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE THRCOR (IC3D, IS, BETA, TCOR, QTOT)
C
C     THRCOR - compute threshold corrections for ICVT
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     EXTRAS - print out extra info about GTS
C     GTST   - compute free energies, CVT and ICVT rates
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     EBEND  - compute bending energy levels
C     ESTR   - compute stretch vibrational energy levels
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LFLG
      PARAMETER (NSTR=10, NBEND=25, NBM1=NBEND-1)
      DIMENSION EVIB(NSTR), EB(NBEND), DEGEN(NBEND)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /COBND/  NLEVEL
      LOGICAL LMAX
      COMMON /EBND1/  LMAX
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/  LGS2
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      TCOR = 1.D0                                                       GCL1092
C     WRITE (99, 9900) IC3D, IS, SS(IS), BETA, QTOT
C9900 FORMAT (' ENTER THRCOR, IC3D,IS,S,BETA,QTOT=', 2I5, 1P3E13.5)
      IF (IS .LT. 1 .OR. IS .GT. NSDM10) RETURN
      VDIF = VMAX(IC3D) - VS(IS)
      IF (VDIF .LE. 0.D0) RETURN
      D = DS(IS)
      XK = XKS(IS)
      E = ESTR(IS, SS(IS), D, XK)
C     WRITE (99, 9901) VS(IS), VMAX(IC3D), VDIF, E
C9901 FORMAT (' VS,VMAX,VDIF,E=', 1P4E13.5)
      IF (E .GT. VDIF) RETURN
      NMAX = 100
      IF (LGS(3)) NMAX = MIN(100, INT((XK-1.D0)/2.D0 + 1.D0))
C     WRITE (99, 9902) NMAX
C9902 FORMAT (' NMAX=', I5)
      IF (NMAX .LT. 1) RETURN
      IF (IC3D .EQ. 1) THEN
C  collinear correction
         T2 = EXP(-BETA*VDIF-QTOT)
         TCOR = TCOR + T2 - EXP(-BETA*E-QTOT)
         I = 1
         IF (LSTATE. EQ. 0 .AND. NMAX .GE. 2) THEN
            TI = 1.D0/XK
            T = 2.D0*D*TI
            TI = 0.5D0*TI
            TIP1 = 1.0D0
            DO 10 I = 2,NMAX
               TIP1 = TIP1 + 2.D0                                       GCL1092
               E = T*TIP1
               IF (LGS(3)) E = E*(1.D0 - TIP1*TI)
               IF (E .GT. VDIF) GO TO 20
               TCOR = TCOR + T2 - EXP(-BETA*E-QTOT)
   10       CONTINUE
   20       CONTINUE
         END IF
         IF (LSTATE .EQ. 0 .AND. E .LE. VDIF) WRITE (6, 6001) I
      ELSE
C  3D correction
         NMAX = MIN(NMAX,NSTR)
         N1 = 1
         EVIB(1) = E
         IF (LSTATE .EQ. 0 .AND. NMAX .GE. 2) THEN
            TI = 1.D0/XK
            T = 2.D0*D*TI
            TI = 0.5D0*TI
            TIP1 = 1.0D0
            DO 30 I = 2,NMAX
               TIP1 = TIP1 + 2.D0                                       GCL1092
               E = T*TIP1
               IF (LGS(3)) E = E*(1.D0 - TIP1*TI)
               IF (E .GE. VDIF) GO TO 40
               N1 = N1 + 1
               EVIB(N1) = E
   30       CONTINUE
   40    CONTINUE
         END IF
C        WRITE (99, 9904) N1, (EVIB(I), I=1,N1)
C9904 FORMAT (' N1=', I5, ', EVIB='/ (1X, 1P10E13.5))
         LFLG = LSTATE .EQ. 0 .AND. E. LT. VDIF
         FB = FBS(IS)
         GB = GBS(IS)
         AB = QFBS(IS)
         IF (LSBEND .EQ. 0) THEN
            IF (LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
               N2 = 0
               EBMX = -1.D10                                            GCL1092
C  Sum over n and K where n = 2v + K
               DO 310 N = 1,NLEVEL
                  IF (N2.GE.NBEND) GO TO 320
                  K = N+2
305               CONTINUE
                     K = K - 2
                     NU = (N-K)/2
                     CALL COBEND(NU, K, FB, AB, GB, EB1)
                     IF (.NOT.LMAX) THEN
                        N2 = N2 + 1
                        DEGEN(N2) = 2.0D0                               GCL1092
                        IF (K.EQ.0) DEGEN(N2) = 1.0D0                   GCL1092
                        EB(N2) = EB1
                        EBMX = MAX(EB1,EBMX)
                     END IF
                  IF (K.GT.1) GO TO 305
310            CONTINUE
320            CONTINUE
               LFLG = LFLG .OR. LMAX .OR. EVIB(1)+EBMX .LT. VDIF
               IF (N2.EQ.0) GO TO 200
               EB2 = EB(1)
            ELSE
C  uncoupled bending energy level
               IF (LGS(4)) THEN
                  EB(1) = EBEND(0, FB, AB, GB)
                  EZ = EB(1) + EVIB(1)
                  E = EZ + EB(1)
C              WRITE (99, 9905) EB(1), EZ, E, VDIF
C9905 FORMAT (' EB(1),EZ,E,VDIF=', 1P4E13.5)
                  IF (E .GE. VDIF) GO TO 200
                  N2 = 1
                  DO 50 I= 1,NBM1
                     E = EBN(I, FB, AB, GB)
                     IF (LMAX .OR. E+EZ .GE. VDIF) GO TO 60
                     N2 = N2 + 1
                     EB(N2)= E
50                CONTINUE
60                CONTINUE
                  LFLG = LFLG .OR. LMAX .OR. E+EZ .LT. VDIF
               ELSE
                  E = VDIF + 1.D0                                       GCL1092
                  IF (FB .LE. 0.0D0) GO TO 200
                  WB = SQRT(FB*GB)
                  EZ = EVIB(1) + 0.5D0*WB
                  N2 = 0
                  E = -0.5D0*WB
                  DO 70 I = 1,NBEND
                     E = E + WB
                     IF (EZ+E .GE. VDIF) GO TO 80
                     N2 = N2 + 1
                     EB(N2) = E
70                CONTINUE
80                CONTINUE
                  E = EZ + E
                  LFLG = LFLG .OR. E .LT. VDIF
                  IF (N2 .EQ. 0) GO TO 200
               END IF
               EB2 = 2.0D0*EB(1)
            END IF
C           WRITE (99, 9906) N2, (EB(I), I=1,N2)
C9906 FORMAT (' N2=', I5, ', EB='/ (1X, 1P10E13.5))
         ELSE
            E = VDIF + 1.D0                                             GCL1092
            IF (LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
               CALL COBEND (NBND1,NBND2,GB,AB,GB,EB2)
            ELSE
C  uncoupled bending energy level
               EB2 = EBEND(NBND1, FB, AB, GB)
               IF (LMAX) GO TO 200
               IF (NBND1 .EQ. NBND2) EB2 = 2.D0*EB2
               IF (NBND1 .NE. NBND2) EB2 = EB2 + EBN(NBND2, FB, AB, GB)
            END IF
            IF (LMAX) GO TO 200
         END IF
         TB = 1.0D0/XMOMS(IS)
         EZ = EVIB(1) + EB2
         T1 = EXP(-BETA*VDIF-QTOT)
         T2 = -1.0D0
         EJ = 0.0D0
         J = -1
         DO 150 IJ = 1,100
            J = J + 1
            EJ = EJ + DBLE(J)*TB                                        GCL1092
            E = EZ + EJ
            IF (E .GE. VDIF) GO TO 160
            T2 = T2 + 2.0D0
            T3 = 0.0D0
            DO 130 I1 = 1,N1
               E1 = EVIB(I1) + EJ
               E = E1 + EB2
               IF (E .GE. VDIF) GO TO 140
               IF (LSBEND .EQ. 0) THEN
                  IF (LGS2(6)) THEN
C  centrifugal oscillator energy levels for bend
                     DO 330 I2 = 1,N2
                        E = E1 + EB(I2)
                        IF (E.LT.VDIF) T3 = T3 + (T1-EXP(-BETA*E-QTOT))
     *                     *DEGEN(I2)
330                  CONTINUE
                  ELSE
C  uncoupled bend energy levels
                     DO 110 I2 = 1,N2
                        EB1 = EB(I2)
                        E2 = EB1 + E1
                        E = E2 + EB1
                        IF (E .GE. VDIF) GO TO 120
                        T3 = T3 + T1 - EXP(-BETA*E-QTOT)
                        IF (I2 .NE. N2) THEN
                           I2P1 = I2 + 1
                           DO 90 I2P = I2P1,N2
                              E = EB(I2P) + E2
                              IF (E .GE. VDIF) GO TO 100
                              T3 = T3 + 2.0D0*(T1 - EXP(-BETA*E-QTOT))
90                         CONTINUE
100                        CONTINUE
                        END IF
110                  CONTINUE
120                  CONTINUE
                  END IF
               ELSE
                  E = E1 + EB2
                  IF (E .GE. VDIF) GO TO 140
                  T3 = T3 + DGBND*(T1 - EXP(-BETA*E-QTOT))
               END IF
  130       CONTINUE
  140       CONTINUE
C           WRITE (99, 9910) J, EJ, T2, T3
C9910 FORMAT (' J,EJ,T2,T3=', I5, 1P3E13.5)
            IF (T3 .EQ. 0.0D0) GO TO 200
            TCOR = TCOR + T2*T3
  150    CONTINUE
  160    CONTINUE
  200    CONTINUE
         LFLG = LFLG .OR. E .LT. VDIF
         IF (LFLG) WRITE (6, 6000) SS(IS), N1, N2, J
      END IF
      IF (TCOR .LE. 0.D0) WRITE (6, 6002) IC3D, IS, SS(IS), TCOR
      RETURN
 6000 FORMAT (2X,T5,'At s=', F10.5,' ICVT correction - (N1,N2,J)=',3I5)
 6001 FORMAT (2X,T5,'ICVT correction - N1(collinear) =', I4)
 6002 FORMAT (2X,T5,'Warning: In THRCOR for IC3D, IS=', 2I4, 
     *       ' and s=', F10.4,', the threshold correction is', 1PE13.5)
      END
