!!!*************************************************************
! 文件/File: rpupdt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpupdt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPUPDT(DELS,DELSS,DELSV2,DELC,DX,DY,DXASY,DYASY,SGN,
     *   SCUR,SRP,VASY,XKASY,VTEST,XKTEST,NSTAB,IC,IEND,IS,
     *   ISGN,ISAD,ISADL,ISMN,ISMX,LASY,LBACK,LEND,LTERM)
C
C     RPUPDT - update S and store rp info
C
C  Called by:
C     RPSEG  - compute reaction path for a segment
C
C  Calls:
C     CURVE  - compute curvature by finite difference
C     RPSTOR - store rp info and check grid
C     SURFIT - calls PARAM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LASY, LBACK, LCUR, LEND, LSAVE, LTERM
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /RPWCOM/ XWL, YWL                                          GCL1096
C
      S = S + DELSS*SGN
      DELSS = 0.D0                                                      GCL1092
C  calculate surface parameters
      CALL SURFIT (DX, DY, SGN)
      LSAVE = .TRUE.
C  decide whether to end this segment
      LEND = (S*SGN-SRP) .GT. -0.0001D0*DELS .AND. .NOT.LBACK .AND.     GCL1092
     * IEND .NE. 0
      IF (.NOT.LEND .AND. IEND .NE. 0 .AND. LASY) THEN
         VDIF = ABS(VNOW-VASY)
         IF (ABS(VASY)  .GT.  1.D-9) VDIF = VDIF/ABS(VASY)              GCL1092
         IF (VDIF .LE. EPSASY) THEN
            XKDIF = ABS((XKM-XKASY)/XKASY)
            IF (XKDIF  .LE.  EPSASY)  THEN
               IF (ABS(FB)  .LT.  EPSASY) LEND = .TRUE.
            END IF
         END IF
C  decide whether to save information and print it out
         IF (.NOT.LEND) THEN
            ISKIP = ISKIP + 1
            IF (ISKIP .LT. 5 .AND. ISADL .GT. 0) THEN
               VDIF = ABS(VOLD-VNOW)
               IF (VDIF .LT. VTEST) THEN
                  XKDIF = ABS(XKOLD-XKM)
                  LSAVE = XKDIF .GT. XKTEST
               END IF
            END IF
         END IF
      END IF
      IF (LSAVE) THEN
         VOLD = VNOW
         XKOLD = XKM
C  check if curvature is to be calculated
         LCUR = S*SGN .LT. SCUR
         IF (LBACK) LCUR = .NOT.LCUR
         LCUR = LCUR .AND. .NOT.LASY .AND. NSTAB.GT.1
         CUR = 0.D0                                                     GCL1092
         IF (LCUR) CALL CURVE (X, Y, DX, DY, DELC, SGN, CUR)
         ISKIP = 0
         RB = Y/CM2
         RA = X - CM1*RB
C  print and store rp info
         WRITE (6, 616) S, RA, RB, X, Y, VNOW, DX, DY, CUR,
     *    LASY, NSTAB
         CALL RPSTOR(ISAD,ISGN,IS,ISMN,ISMX,IEND,LTERM)                 GCL0493
      END IF
      IF (LEND .OR. LTERM) RETURN
C  reset grid save size value and step size for curvature
      IF (ISADL .LT. 3) THEN
         DELSV2 = DELSV*.399999D0                                       GCL1092
         ISADL = ISADL + 1
         IF (ISADL .GT. 2) THEN
            DELSV2 = DELSV*0.999999D0                                   GCL1092
            DELC = DELCUR
         END IF
      END IF
C  reset NSTAB and possibly LASY
      NSTAB = MAX(NSTAB, 10)
      IF (IEND .EQ. 0) THEN
         LASY = .FALSE.
      ELSE IF (LBACK) THEN
C  For LBACK = .TRUE. reset LASY to false and check the
C     distance to well
         DXWL = X - XWL
         DYWL = Y - YWL
         DMAG = SQRT(DXWL*DXWL + DYWL*DYWL)
         LASY = .FALSE.
         IF (DMAG .LE. DELSV2) NSTAB = 1
      ELSE IF (.NOT. LASY .AND. S*SGN .GT. SCUR) THEN
C  For LBACK = .FALSE., if LASY = .FALSE. and S is outside the region
C     where the curvature is computed, test if LASY can be set true.
         AC = DX*DXASY + DY*DYASY
         IF (AC .GE. ACASY) THEN
            IC = IC + 1
            IF (IC .GT. 3) THEN
               DX = DXASY
               DY = DYASY
               LASY = .TRUE.
               NSTAB = 10
            END IF
         ELSE
            IC = 0
         END IF
      END IF
      RETURN
C
  616 FORMAT(1X,F10.6,2(2X,F10.6,1X,F10.6),1X,1P,E13.5,2X,0P,E12.5,1X,  TCA1097
     &       E12.5,1X,1P,E12.4,3X,0P,L1,3X,I5)                          TCA1097
      END
