!!!*************************************************************
! 文件/File: rpseg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpseg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPSEG (DELFIX,ERREPS,VTEST,XKTEST,ISWL,IEND,ISAD,ISGN,
     *                  IS,ISMN,ISMX,LTERM,LRFLC)
C
C     RPSEG  - compute reaction path for a segment
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     INIT   - initialization of arrays for following reaction path
C     RPSEG2 - compute rp to next save point
C     RPSTOR - store rp info and check grid
C     RPUPDT - update S and store rp info
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LASY, LBACK, LEND, LTERM, LWELL, LRFLC                    GCL1096
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      PARAMETER (NSADDM=4)
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     *OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /RP2/    SLM, SLP
      COMMON /RPWCOM/ XWL, YWL                                          GCL1096
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)

C  initialization to saddle point values
      CALL INIT (ISAD)
      SGN = DBLE(ISGN)                                                  GCL1092
      DX = -SGN*SVECT(1,ISAD)
      DY = -SGN*SVECT(2,ISAD)
      WRITE (6, 614) ISAD
      WRITE (6, 616) S, R1SP(ISAD), R2SP(ISAD), X, Y, VNOW, DX, DY
C  store saddle point info
      NSSP(ISAD) = IS
      CALL RPSTOR(ISAD,ISGN,IS,ISMN,ISMX,IEND,LTERM)                    GCL0493
      IF (LTERM) RETURN
C  reset variables
      IC = 0
      ISADL = 0
      NSTAB = 10
      LASY = .FALSE.
      LBACK = .FALSE.
      LEND = .FALSE.
      LWELL = .FALSE.
      DELC = DELCSP
      DELSS = 0.D0                                                      GCL1092
      DELSV2 = DELSV*0.199999D0                                         GCL1092
C  for IEND .NE. 0, set up asymptotic values
      IF (IEND .EQ. 1) THEN
C  products
         VASY = VDELTA
         XKASY = XKAB
         T = SQRT(CM1*CM1+CM2*CM2)
         DXASY = -CM1/T
         DYASY = -CM2/T
         SCUR = SP
         SRP = SLP
      ELSE IF (IEND .EQ. -1) THEN
C  reactants
         VASY = 0.D0                                                    GCL1092
         XKASY = XKBC
         DXASY = -1.0D0                                                 GCL1092
         DYASY = 0.D0                                                   GCL1092
         SCUR = -SM
         SRP = -SLM
      END IF
C  loop over save points
10    CONTINUE
C  reset step size and save size
         DEL = DELFIX
         DELSV3 = DELSV2
C  compute rp to next save point
         CALL RPSEG2 (DELS,DELSS,DELSV2,DELSV3,DX,DY,DXASY,DYASY,
     *                ERREPS,SGN,VASY,SRP,ISWL,NSTAB,IEND,IS,ISGN,      GCL0493
     *                ISAD,ISADL,LASY,LBACK,LEND,LRFLC,LTERM,LWELL)     GCL1096
         IF (.NOT. (LEND .OR.LTERM))
C  update S and store rp info
     *      CALL RPUPDT(DELS,DELSS,DELSV2,DELC,DX,DY,DXASY,DYASY,SGN,
     *        SCUR,SRP,VASY,XKASY,VTEST,XKTEST,NSTAB,IC,IEND,           GCL1096
     *        IS,ISGN,ISAD,ISADL,ISMN,ISMX,LASY,LBACK,LEND,LTERM)
      IF (.NOT. (LEND .OR.LTERM)) GO TO 10
      RETURN
  614 FORMAT (1X,T2,'Saddle point', I3)
  616 FORMAT(1X,F10.6,2(2X,F10.6,1X,F10.6),1X,1P,E13.5,2X,0P,E12.5,1X,  TCA1097
     &       E12.5)                                                     TCA1097
      END
