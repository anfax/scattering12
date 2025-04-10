!!!*************************************************************
! 文件/File: rpwell.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpwell.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPWELL(NSTAB,DELSS,DELSV2,DX,DY,DXASY,DYASY,XWLSV,
     *                  YWLSV,SGN,VASY,SRP,ISWL,IEND,IS,ISGN,ISAD,
     *                  LASY,LBACK,LEND,LRFLC,LTERM) 
C
C     RPWELL - check distance to well
C
C  Called by:
C     RPSEG2 - compute rp to next save point
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LASY,LBACK,LEND,LRFLC,LTERM
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /RPWCOM/ XWL, YWL                                          GCL0693
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      DXWL = X - XWL
      DYWL = Y - YWL
      DMAG = SQRT (DXWL*DXWL + DYWL*DYWL)
      IF (DMAG .LE. DEL) THEN
C  less than DEL from minimum, take last step
         DELSS = DELSS + DMAG
         IF (DELSS .LT. .5D0*DELSV2) IS = IS - ISGN                     GCL1092
         X = XWL
         Y = YWL
C  at bottom of well, print out rp info about well
         RB = Y/CM2
         RA = X - CM1*RB
         S = S + SGN*DELSS
         DELSS = 0.D0                                                   GCL1092
         CUR = 0.D0                                                     GCL1092
         WRITE (6, 616) S, RA, RB, X, Y, VNOW, DX, DY, CUR, LASY, NSTAB
         IF ((IEND .NE. 0 .OR. ISGN .GT. 0) .AND. .NOT.LBACK) THEN
C  at well first time, write out special information about the well,
C     then store rp info
            IF (ABS(S-SS(IS-ISGN)) .LT. 0.5D0*DELSV2) IS = IS - ISGN    GCL1092
            WRITE (6, 618)
            CALL WELL (XWL, YWL, UXWL, UYWL, IS)
            DX = SGN*UYWL
            DY = -SGN*UXWL
            X = XWL
            Y = YWL
            XWLSV = XWL
            YWLSV = YWL
            CALL SURFIT (DX, DY, SGN)
            CALL STORE (IS)
            ISWL = IS
C  if symmetric system, then done with this segment
            IF (LSYM .AND. ISGN .GT. 0 .AND. ABS(RA-RB) .LT. 1.D-8) THEN GCL1092
               IS = IS + ISGN
               LRFLC = .TRUE.
               LEND = .TRUE.
               RETURN
C  if IEND = 0, then done with this segment
            ELSE IF (IEND .EQ. 0) THEN
               IS = IS + ISGN
               LEND = .TRUE.
               RETURN
            END IF
C  make sure the well is below the asymptote, otherwise another saddle
C     point exits
            IF (VNOW .GE. VASY) THEN
               WRITE (6, 6005)
               STOP 'RPWELL 1'
            END IF
C  set IS and ISMN or ISMX and check if theres's enough room left in
C     grid for more points
            IF (ISGN .GT. 0) THEN
               ISMN = ISWL +1
               IF (ISMN .GT. NSMAX) THEN
                  LTERM = .TRUE.
                  RETURN
               END IF
               IS = NSMAX
            ELSE
               ISMX = ISWL - 1
               IF (ISMX .LT. 1) THEN
                  LTERM = .TRUE.
                  RETURN
               END IF
               IS = ISMN
            END IF
C  set up to start in asymptotic region and integrate towards well
            DXASY = -DXASY
            DYASY = -DYASY
            DX = DXASY
            DY = DYASY
            DXSV = DX
            DYSV = DY
C  move out into asymptotic well
10          CONTINUE
               IF (ISGN .LT. 0) RA = RA + 2.D0                          GCL1092
               IF (ISGN .GT. 0) RB = RB + 2.D0                          GCL1092
               X = RA + CM1*RB
               Y = CM2*RB
C  stabilize for new point
               CALL VMIN (X, Y, DX, DY, DEL, VNOW, .FALSE., IERR)
               IF (IERR .NE. 0) THEN
                  WRITE (6, 6006)
                  STOP 'RPWELL 2'
               END IF
               T1 = X - XWL
               T2 = Y - YWL
               SNEW = S + SGN*SQRT(T1*T1+T2*T2)
               T = VNOW-VASY
               IF (ABS(VASY)  .GT.  EPSASY) T = T/VASY
            IF (SNEW*SGN.LE.SRP .AND. ABS(T).GT.EPSASY) GO TO 10
C  reset variables for the segment from the asymptote back to well and
C     store rp info for the current location
            ISGN = -ISGN
            NSTAB = 10
            LASY = .TRUE.
            LBACK = .TRUE.
            CUR = 0.D0                                                  GCL1092
            DELSV2 = DELSV*0.999999D0                                   GCL1092
            S = SNEW
            SGN = -SGN
            WRITE (6, 616) S, RA, RB, X, Y, VNOW, DX, DY, CUR,
     *       LASY, NSTAB
            CALL SURFIT (DX, DY, SGN)
            CALL RPSTOR(ISAD,ISGN,IS,ISMN,ISMX,IEND,LTERM)              GCL0493
         ELSE IF (IEND .NE. 0) THEN
C  at well second time, shift grid values so the two segments are
C     contiguous in the grid
            SDIF = SS(ISWL) - S
            WRITE (6, 620) SDIF
            IF (ABS(S-SS(IS-ISGN)) .LT. 0.5D0*DELSV2) IS = IS - ISGN    GCL1092
            IS = IS - ISGN
            ISHIFT = ISWL - IS - ISGN
            CALL SHIFT (IS, ISHIFT, ISGN, SDIF, ISAD, NSSP)
C  reset IS, ISGN and end this segment
            IF (ISGN .GT. 0) THEN
               IS = ISWL - IS - 1
            ELSE
               IS = ISWL + NSMAX - IS + 2
            END IF
            SGN = -SGN
            ISGN = -ISGN
            LEND = .TRUE.
         ELSE
            LEND = .TRUE.
         END IF
      ELSE IF (DMAG .LE. DELSV2) THEN
C  close to bottom of well, go straight towards it
         DX = DXWL/DMAG
         DY = DYWL/DMAG
      END IF
  616 FORMAT(1X,F10.6,2(2X,F10.6,1X,F10.6),1X,1P,E13.5,2X,0P,E12.5,1X,  TCA1097
     &       E12.5,1X,1P,E12.4,3X,0P,L1,3X,I5)                          TCA1097
  618 FORMAT (1X,T5,'At bottom of well')
  620 FORMAT (1X,T5,'Previous s values are shifted by', F8.4)
 6005 FORMAT (/,1X,T5,'Error: The minimum of the well is greater than ',
     * 'the asymptotic value, therefore there is another saddle point.')
 6006 FORMAT (/,1X,T5,'Error: In RPATH there is a problem with VMIN ',
     * 'in finding asymptotic value of reaction path for LROT=.FALSE.')
      END
