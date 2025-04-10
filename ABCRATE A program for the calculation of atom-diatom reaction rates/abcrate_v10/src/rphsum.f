!!!*************************************************************
! 文件/File: rphsum.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rphsum.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPHSUM
C
C     RPHSUM - truncates and summarize reaction path info
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     CURVEX - extrapolate curvature
C     EBEND  - compute bending energy levels
C     RPWRIT - write out reaction path info
C     SHIFT  - shift r.p. info in grid
C     WKBSET - set up grid of WKB energy levels
C     PEF    - evaluate potential
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LPRNT
      DIMENSION VT(9)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /RP2/    SLM, SLP
      COMMON /RP3/    NPRP
      COMMON /RP4/ STRUNL, STRUNR
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      WRITE (6, 600)
C
C  write out information to save file if not restarting
      IF (LGS2(7) .AND. .NOT.LRSTRT) CALL RPWRIT
      N = NSADMX(2)
      IS = NSSP(N)
      IPRNT = NPRP
      ISAD = 0
      ISP = 0
      SDIF = -SS(IS)
      WRITE (6, 602) SDIF
      IF (ABS(SDIF) .GT. 1.D-8) THEN
         DO 10 IS = 1,NSMAX
            SS(IS) = SS(IS)+SDIF
   10    CONTINUE
      END IF
      IF (STRUNL .EQ. 0.0D0 .AND. STRUNR .EQ. 0.0D0) THEN               GCL1092
         STRUNL = SLM
         STRUNR = SLP
      END IF
      NSADO = NSAD
      ST = MAX(SLM,STRUNL) - 1.D-4
      IS = 1
   20 CONTINUE
         IF (IS .GT. NSMAX .OR. SS(IS) .GT. ST) GO TO 30
         IS = IS + 1
      GO TO 20
   30 CONTINUE
      IS = IS - 1
      IF (IS .GT. 0) THEN
         IF (IS .GE. NSMAX) THEN
            WRITE (6, 6000)
            STOP 'RPHSUM 1'
         END IF
         ISHIFT = -IS
         CALL SHIFT (NSMAX, ISHIFT, 1, 0.0D0, NSAD, NSSP)
         NSMAX = NSMAX - IS
         NMID = NMID - IS
         ISAD = NSAD
   40    CONTINUE
            IF (ISAD .LE. 0 .OR. NSSP(ISAD) .LE. 0) GO TO 50
            ISAD = ISAD - 1
         GO TO 40
   50    CONTINUE
         IF (ISAD .GT. 0) THEN
            NSAD = NSAD - ISAD
            I = 1
   60       CONTINUE
               IF (I .GT. NSAD) GO TO 70
               NSSP(I) = NSSP(I + ISAD)
               I = I + 1
            GO TO 60
   70       CONTINUE
         END IF
      END IF
      ST = MIN(SLP,STRUNR) + 1.D-4
   80 CONTINUE
         IF (NSMAX .LE. 0 .OR. SS(NSMAX) .LT. ST) GO TO 90
         NSMAX = NSMAX - 1
      GO TO 80
   90 CONTINUE
      IF (NSMAX .LE. 0) THEN
         WRITE (6, 6000)
         STOP 'RPHSUM 2'
      END IF
  100 CONTINUE
         IF (NSAD .LE. 0 .OR. NSSP(NSAD) .LE. NSMAX) GO TO 110
         NSAD = NSAD - 1
      GO TO 100
  110 CONTINUE
      IF (NSAD .NE. NSADO) THEN
         IF (NSAD .LE. 0) THEN
            NSADMX(1) = 1
            NSADMX(2) = 1
            IF (SLP .LT. 0.0D0) THEN
               NSSP(1) = 1
               NSSP(2) = 1
            ELSE IF (SLM .GT. 0.0D0) THEN
               NSSP(1) = NSMAX
               NSSP(2) = NSMAX
            END IF
         ELSE
            IF (LGS(7)) THEN
               DO 120 I = 1,NSAD
                  IS = NSSP(I)
                  ISS = -2-I
                  E = DS(IS)*TNP1*(2.0D0-TNP1/XKS(IS))/XKS(IS)
                  CALL WKBSET (ISS, SS(IS), NSTATE, E, UL, UG, PER)
                  IF (LSTATE .NE. 0 .AND. NSTATE. NE. 0) THEN
                     E = DS(IS)*(2.0D0-1.0D0/XKS(IS))/XKS(IS)
                     CALL WKBSET (ISS, SS(IS), 0, E, UL, UG, PER)
                  END IF
  120          CONTINUE
            END IF
            DO 140 IC3D = 1,2
               NSADMX(IC3D) = NSADMX(IC3D) - ISAD
               IF (NSADMX(IC3D) .LE. 0) THEN
                  IF (NSAD .EQ. 1) THEN
                     NSADMX(IC3D) = 1
                  ELSE
                     VMX = -1.D30
                     NMX = 1
                     DO 130 I = 1,NSAD
                        IS = NSSP(I)
                        ISS = -2-I
                        VAD = VSP(IS) + 
     *                        ESTR(ISS, SS(IS), DS(IS), XKS(IS))        GCL0795
                        IF (IC3D .EQ. 2) THEN
                           IF (LGS2(6)) THEN
C  centrifugal oscillator energy level
                              CALL COBEND(0,0,FBS(IS),QFBS(IS),
     *                           GBS(IS),EB)
                           ELSE
C  uncoupled bending energy level
                              EB = 2.0D0*EBEND(0,FBS(IS),QFBS(IS),      GCL1092
     *                                         GBS(IS))                 GCL1092
                           END IF
                        END IF
                        IF (VAD .GE.  VMX .AND.  ABS(VAD-VMX) .GE.
     *                     1.D-6) THEN
                           NMX = I
                           VMX = VAD
                        END IF
  130                CONTINUE
                     NSADMX(IC3D) = NMX
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
C  Write out summary information
      IF (LSYM) THEN
         NSLAST = NMID
      ELSE
         NSLAST = NSMAX
      END IF
      DO 150 IS = 1,NSLAST
         IF (IS .GE. ISP) THEN
            ISAD = ISAD + 1
            IF (ISAD .LE. NSAD) THEN
               ISP = NSSP(ISAD)
            ELSE
               ISP = NSLAST
            END IF
            LPRNT = .TRUE.
         ELSE
            IPRNT = IPRNT + 1
            LPRNT = IPRNT .GE. NPRP
         END IF
         IF (LPRNT) THEN
            IPRNT = 0
            V = VS(IS)*CKCAL
            D = DS(IS)
            XK = XKS(IS)
            IF (XK .NE. 0.D0) THEN
               WSTR = CCM * 4.D0 * D / XK
            ELSE
               WSTR = 1.D35
            END IF
            D = D * CKCAL
            FB = FBS(IS)
            QFB= QFBS(IS) * CKCAL
            GB = GBS(IS)
            WB = 0.D0                                                   GCL1092
            IF (FB .GT. 0.D0) WB = CCM*SQRT(GB*FB)
            FB = FB*CKCAL
            WRITE (6, 604) SS(IS), V, D, XK, WSTR, FB, QFB, GB, WB,
     *         XMOMS(IS)
         END IF
  150 CONTINUE
      IF (LSYM) WRITE (6, 610)
      IF (.NOT.LRSTRT) THEN
C  test of bending potential
         WRITE (6, 606)
         DELPHI = PI*5.D0/180.D0
         DO 170 IS = 1,NSLAST,10
            X = X2(IS)
            Y = Y2(IS)
            R2 = Y/CM2
            R1 = X - CM1*R2
            PHI = PI - DELPHI
            DO 160 I = 1,9
               PHI = PHI + DELPHI
               R3 = SQRT(R1*R1 + R2*(R2-2.D0*R1*COS(PHI)))
               CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                  GCL0893
               VT(I) = V*CKCAL
  160       CONTINUE
            WRITE (6, 608) SS(IS), R1, R2, VT
  170    CONTINUE
      END IF
      CALL CURVEX
      RETURN
  600 FORMAT (//, T36, 'Morse parameters', T78, 'Bending parameters',
     *   T115, 'Moment',/, 4X, 's', T17, 'V', T31, 'De', T42, 'K', T48,
     *   'hbar W(STR)', T64, 'FPHI', T77, 'APHI', T90, 'GPHI', T99,
     *   'hbar W(BND)', T113, 'of inertia', /, 2X,                      TCA1097
     *   '(bohr)', T13, '(kcal)', T27, '(kcal)', T49,
     *   '(cm**-1)', T61, '(kcal)', T74, '(kcal)', T89,
     *   '(a.u.)', T100, '(cm**-1)', T115, '(a.u.)')
  602 FORMAT (1X,T5,'Note: Previous s values are shifted by', F8.4)
  604 FORMAT (1X, F8.4, 1P,2E14.5, 0PF9.4, F10.2, 2X, 1P,3E13.5,  
     1   0PF10.2, 1PE16.5)                                              TCA1097
  606 FORMAT (/,1X,T5,'Test of the bending potential',//,72X,
     *   'V(kcal)'/ T6, 's', T14, 'R1', T21, 'R2', T26, 'PHI=', T30,
     *   '0', T41, '5', T51, '10', T63, '15', T76, '20', T86, '25', T96,
     *   '30', T107, '35', T118, '40'/)
  608 FORMAT (1X, F7.3, 1X, 2F7.3, 1X, 1P,9E11.3)     
  610 FORMAT (/,2X,T5,'Note: Symmetric system, product side ',
     *                'obtained by reflection')
 6000 FORMAT (/,2X,T5,'Error: The grid has been shifted such that no ',
     *                'entries remain.')
      END
