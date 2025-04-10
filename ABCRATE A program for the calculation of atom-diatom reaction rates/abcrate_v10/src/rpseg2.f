!!!*************************************************************
! 文件/File: rpseg2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpseg2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPSEG2 (DELS,DELSS,DELSV2,DELSV3,DX,DY,DXASY,DYASY,
     *                   ERREPS,SGN,VASY,SRP,ISWL,NSTAB,IEND,IS,ISGN,
     *                   ISAD,ISADL,LASY,LBACK,LEND,LRFLC,LTERM,LWELL)
C
C     RPSEG2 - compute rp to next save point
C
C  Called by:
C     RPSEG  - compute reaction path for a segment
C
C  Calls:
C     GRAD   - follow gradient
C     INIT   - initialization of arrays for following reaction path
C     ROOT2D - 2d root search
C     RPWELL - check distance to well
C     STEP1  - take step along a vector then find minimum normal to
C              vector
C     STEP2  - follow gradient by choosing step to minimize potential
C              along an arc
C     VMIN   - find minimum energy along u coordinate
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DERIV
      LOGICAL LASY,LBACK,LEND,LRFLC,LTERM,LWELL
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     * CUR
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     * NSSP(NSADDM), LRSTRT
      COMMON /RPWCOM/ XWL, YWL                                          GCL0693
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     * XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     * NSAD, NSADMX(2)
C
      LEND = .FALSE.
C  loop over NSTAB steps (or one step of NSTAB*DEL)
10    CONTINUE
C        WRITE (91, 9103) LASY, NSTAB
C9103    FORMAT (' LASY,NSTAB=', 1X, L1, I5)
      DELN = DEL * NSTAB
      IF (LASY) THEN
C  LASY = .TRUE. (in asymptotic region), take one step of NSTAB*DEL
         STEP = MIN(DELN, DELSV-DELSS)
         STEP = MAX(STEP, DEL)
         CALL STEP1 (X, Y, DX, DY, STEP, VNOW, DELS, ERR)
         DELSS = DELSS + DELS
         IF (ERR .GE. DEL) NSTAB = MAX(1, NSTAB-10)
         IF (ERR .LT. ERREPS) NSTAB = NSTAB+5
         DELN = DEL*NSTAB
      ELSE IF (NSTAB .GT. 1) THEN
C  follow gradient, repeatedly take NSTAB steps of DEL until DELSS >
C     DELSV2
         CALL GRAD (X, Y, DX, DY, DEL, DELS, DELSS, DELSV2, VNOW,
     *    NSTAB, IRET)
C  IRET = 1, normal completion
C  IRET = 2, too many oscillations
         IF (IRET .EQ. 1) GO TO 999
         NSTAB = 1
         IF (IEND .NE. 0) THEN
C  check if in asymptotic region
            AC = DX*DXASY + DY*DYASY
            IF (AC .GE. ACASY) THEN
C  assumed in asymptotic region
               LASY = .TRUE.
               DX = DXASY
               DY = DYASY
            END IF
         END IF
C           WRITE (91, 9104) LASY
C9104       FORMAT (' IRET=2, LASY=', 1X, L1)
         IF (.NOT.LASY) THEN
C  check for well
            IF (.NOT.LWELL) THEN
               NMAX = 25
               XWL = X
               YWL = Y
               CALL ROOT2D (DERIV, XWL, YWL, DEL, 1.D-10, 0.0D0,        GCL1092
     *          NMAX)
               LWELL = .FALSE.
               IF (NMAX .EQ. 0) GO TO 10
C  check that the extrema found is not a saddle point
               IF (NSAD .GT. 0) THEN
                  DO 20 I = 1,NSAD
                     IF (ABS(XWL-XSP(I)) .GT. 1.D-6) GO TO 20           GCL1092
                     IF (ABS(YWL-YSP(I)) .LT. 1.D-6) GO TO 10           GCL1092
20                CONTINUE
               END IF
               LWELL = .TRUE.
C                 WRITE (91, 9105)
C9105             FORMAT (' WELL FOUND')
               IF (IEND .EQ. 0 .AND. ISGN .LT. 0)  THEN
C  For IEND=0, ISGN<0, this should be the second time into this well.
C     Check that the new and old coordinates of the well agree.
                  IF (ABS(XWL-XWLSV) .GT. 1.D-6 .OR.                    GCL1092
     *             ABS(YWL-YWLSV) .GT. 1.D-6) THEN                      GCL1092
                     WRITE (6, 6001)
                     STOP 'RPSEG2 1'
                  END IF
               END IF
            END IF
C     Check distance to well)
            CALL RPWELL(NSTAB,DELSS,DELSV2,DX,DY,DXASY,DYASY,XWLSV,
     *                  YWLSV,SGN,VASY,SRP,ISWL,IEND,IS,ISGN,ISAD,      GCL1096
     *                  LASY,LBACK,LEND,LRFLC,LTERM)                    GCL1096
         END IF
      ELSE
C  NSTAB = 1, take steps to lowest energy point along arcs of length
C     STEP
30       CONTINUE
         STEP = MIN(DELN, DELSV3-DELSS)
         STEP = MAX(STEP, DEL)
C              WRITE (91, 9106) STEP, DELSS, DELSV3
C9106          FORMAT (' STEP,DELSS,DELSV3=', 1P3E13.5)
         XSV = X
         YSV = Y
         CALL STEP2 (X, Y, DX, DY, STEP, VNOW)
         T1 = X-XSV
         T2 = Y-YSV
         DELSS = DELSS + SQRT(T1*T1 + T2*T2)
         IF (DELSS .LT. DELSV3) THEN
            IF (LWELL) THEN
C     Check distance to well)
               CALL RPWELL(NSTAB,DELSS,DELSV2,DX,DY,DXASY,DYASY,XWLSV,
     *                     YWLSV,SGN,VASY,SRP,ISWL,IEND,IS,ISGN,ISAD,   GCL1096
     *                     LASY,LBACK,LEND,LRFLC,LTERM)                 GCL1096
               GO TO 99
            END IF
         ELSE
C  ready to store rp info, try to stabilize, first try with LROT =
C  .true.
            XSV = X
            YSV = Y
            DXSV = DX
            DYSV = DY
            V = VNOW
            CALL VMIN (X, Y, DX, DY, STEP, VNOW, .TRUE., IERR)
C                    WRITE (91, 9107) IERR
C9107                FORMAT (' LROT=.T., IERR=',I4)
            IF (IERR .NE. 0) THEN
C  if first try didn't work try again with LROT = .false.
               X = XSV
               Y = YSV
               DX = DXSV
               DY = DYSV
               CALL VMIN (X, Y, DX, DY, STEP, VNOW, .FALSE.,
     *          IERR)
C                       WRITE (91, 9108) IERR
C9108                   FORMAT (' LROT=.F., IERR=',I4)
               IF (IERR .NE. 0) THEN
C  unable to stabilize, write error message
                  WRITE (6, 6002) XSV, YSV, DXSV, DYSV
                  IF (ISADL .EQ. 0) THEN
C  the first step from the saddle point failed, increase step size and
C     try again
C                             WRITE (91, 9109) DELSS, DEL
C9109                         FORMAT (' ISADL=0, DELSS, DEL=', 1P2E13.5)
                     DEL = DEL * 1.5D0                                  GCL1092
                     IF (DEL .GE. DELSV3) THEN
                        WRITE (6, 6003)
                        STOP 'RPSEG2 2'
                     END IF
                     DELSS = 0.D0                                       GCL1092
                     CALL INIT(ISAD)
                     DX = -SGN*SVECT(1,ISAD)
                     DY = -SGN*SVECT(2,ISAD)
                  ELSE
C  reset DELSV3 and try taking some more steps
                     X = XSV
                     Y = YSV
                     DX = DXSV
                     DY = DYSV
                     VNOW = V
                     DELSV3 = DELSV3 + DELSV2
C  if DELSV3 gets too big quit
                     IF (DELSV3 .GT. 20.D0*DELSV2) THEN                 GCL1092
                        WRITE (6, 6004)
                        STOP 'RPSEG2 3'
                     END IF
                  END IF
               END IF
            END IF
         END IF
         IF (DELSS.LT.DELSV3 .AND. .NOT.(LEND.OR.LTERM)) GO TO 30
99       CONTINUE
      END IF
      IF (DELSS.LT.DELSV2 .AND. .NOT.(LEND.OR.LTERM)) GO TO 10
999   CONTINUE
      RETURN
6001  FORMAT (/,1X,T5,'Error: In RPSEG2 for a multiple saddle point ',
     *                'surface the coordinates of the ',
     *        /,1X,T12,'well supposedly connecting two saddle points',
     *                 'do not have the same ',
     *        /,1X,T12,'location coming from the different saddle ',
     *                 'points')
6002  FORMAT (/,1X,T5,'Error: In RPSEG2 problem with VMIN, NSTAB=1, ',
     *                'LROT=.FALSE.',
     *        /,1X,T5,'X,Y,DX,DY=', 1P,4E13.5) 
6003  FORMAT (/,1X,T5,'Error: In RPSEG2 unable to take a large enough ',
     *                'step from saddle point to start following grad.')
6004  FORMAT (/,1X,T5,'Unable to stabilize before storing RP info.',
     * ' quiting after 20 attempts.')
      END 
