!!!*************************************************************
! 文件/File: wkb.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: wkb.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WKB (IS, S, NLVL, EIG, UL, UG, PER)
C
C     WKB    - compute WKB energy levels for stretch
C
C  Called by:
C     LAG    - compute LAG probabilities
C     VSPLN2 - spline fit of adiabatic potential in product channel
C     WKBSET - set up grid of WKB energy levels
C
C  Calls:
C     VEXTR  - find extrema in potential in u coordinate (for WKB
C              quantization)
C     PHSINT - compute phase integrals needed for WKB quantization
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCONV
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /WKB1/   X0, Y0, UX, UY
      COMMON /WKB2/   UMIN, VMIN, UMAX, VMAX
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE SOLD, NLVLO                                                  TCA1097
      DATA SOLD /1.D36/, NLVLO /1000000/
C
C     WRITE (99, 9902) IS, S, SOLD, NLVL, NLVLO
C9902 FORMAT (' WKB: IS,S,SOLD,NLVL,NLVLO=', I5, 1P2E13.5, 2I5)
      IF (S .NE. SOLD .OR. NLVL .NE. NLVLO) THEN
         IF (S .NE. SOLD) THEN
            IF (IS .GT. 0 .AND. IS .LE. NSMAX) THEN
               X0 = X2(IS)
               Y0 = Y2(IS)
               UX = UXS(IS)
               UY = UYS(IS)
            ELSE IF (IS .EQ. -1) THEN
               X0 = 1000.D0                                             GCL1092
               Y0 = CM2*R2ASY
               UX = 0.D0                                                GCL1092
               UY = 1.D0                                                GCL1092
            ELSE IF (IS .EQ. -2) THEN
               Y0 =1000.D0                                              GCL1092
               X0 = R1ASY + CM1*Y0/CM2
               T = SQRT(CM1*CM1+CM2*CM2)
               UX = CM2/T
               UY =-CM1/T
            ELSE IF (IS .LT. -2) THEN
               ISAD = -IS-2                                             GCL1096
               X0 = XSP(ISAD)                                           GCL1096
               Y0 = YSP(ISAD)                                           GCL1096
               UX = UVECT(1,ISAD)                                       GCL1096
               UY = UVECT(2,ISAD)                                       GCL1096
            ELSE
               WRITE (6, 6000) IS
               STOP
            END IF
C           WRITE (99 ,9903) X0, Y0, UX, UY
C9903       FORMAT (' X0,Y0,UX,UY=', 1P4E13.5)
            CALL VEXTR
         END IF
C        WRITE (99, 9904) UMIN, VMIN, UMAX, VMAX
C9904    FORMAT (' UMIN,VMIN,UMAX,VMAX=',1P4E13.5)
         XN = DBLE(NLVL) + 0.5D0                                        GCL1092
         U1 = 0.D0                                                      GCL1092
         U2 = 0.D0                                                      GCL1092
         EMIN = VMIN
         EMAX = VMAX - 0.00000001D0*ABS(VMAX)
         E = EMAX
         CALL PHSINT (E, U1, U2, THETA, DNDE, IERR)
         IF (IERR .NE. 0) WRITE (6, 6003) IS, S, NLVL, E*CKCAL
         IF (THETA .LT. XN) THEN
            WRITE (6, 6001) IS, S, NLVL, VMAX, THETA
         ELSE
            E = EIG + VMIN
            IF (E .LE. EMIN .OR. E .GE. EMAX) E = 0.5D0*(EMIN+EMAX)
            IC = 0
C           WRITE (99, 9900) EMIN, EMAX, THETA, U1, U2, E
C9900       FORMAT (' EMIN,EMAX,THETA(EMAX),U1,U2,E=', 1P6E13.5)
   10       CONTINUE
               IC = IC + 1
               CALL PHSINT (E, U1, U2, THETA, DNDE, IERR)
               IF (IERR .NE. 0) WRITE (6, 6003) IS, S, NLVL, E*CKCAL
               F = THETA - XN
C              WRITE (99, 9901) IC, E, F, U1,U2,THETA, EMIN, EMAX
C9901          FORMAT (' IC,E,F,U1,U2,THETA,EMIN,EMAX=', I5, 1P7E15.7)
               LCONV = ABS(F) .LT. 1.D-8
               IF (.NOT.LCONV) THEN
                  IF (F .LT. 0.0D0) THEN
                     EMIN = E
                  ELSE
                     EMAX = E
                  END IF
                  IF (DNDE .NE. 0.0D0) THEN
                     E = E - F/DNDE
                  ELSE
                     E = EMIN
                  END IF
                  IF (E .LE. EMIN .OR. E .GE. EMAX) E =
     *               0.5D0*(EMIN+EMAX)
                  LCONV = ABS(EMAX-EMIN) .LT. 1.D-9
               END IF
            IF (.NOT.LCONV .AND. IC .LE. 100) GO TO 10
            IF (.NOT.LCONV) WRITE (6, 6004) IS, S, NLVL, E, THETA
         END IF
      END IF
      EIG = E - VMIN
      UL = U1
      UG = U2
      PER = TPI*DNDE
      SOLD = S
      NLVLO = NLVL
      RETURN
 6000 FORMAT (' WKB WAS CALLED WITH IS=', I5)
 6001 FORMAT (' WKB, IS,S=', I5, F15.5, ', FOR STATE', I5,
     *   ', VMAX=', 1PE17.9, ', AT VMAX, THETA=', E17.9)
 6002 FORMAT (' WKB, IS,S=', I5, F15.5, ', FOR STATE', I5,
     *   ', VMIN=', 1PE17.9, ', AT VMIN, THETA=', E17.9)
 6003 FORMAT (' ERROR IN PHSINT FOR IS,S,NLVL,E=', 2(I5, 1PE17.9))
 6004 FORMAT (' ITERATION TO WKB EIGENVALUE NOT CONVERGED, IS,S=', I5,
     *   F15.5, ', FOR STATE', I5, ', E=', 1PE17.9, ', THETA=', E17.9)
      END
