!!!*************************************************************
! 文件/File: lbetas.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lbetas.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LBETAS (E,ALF,BETA,XAB,YAB,X,Y,UX,UY,UALF,DELSMN,
     *   DELSMX,S,ISINT,IERR)
C
C     LBETAS - solve for s of beta
C
C  Called by:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C
C  Calls:
C     LALF0  - set up zeroth order path for LAG
C     LBUANG  - compute UALF and angle between normalized gradient
C               vector at s and vector from xmep at s to point
C               on the tunneling path
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LCONV
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /THLAG1/ XM1, YM1, X1, Y1, UXM1, UYM1, UX1, UY1, SL, SR,
     *XDIF, XSUM, YDIF, YSUM
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE EPSA, EPSD, EPSU, CKCAL                                      TCA1097
      DATA EPSA/1.0D-8/, EPSD/1.0D-8/, EPSU/1.D-9/                      GCL1092
      DATA CKCAL/627.5095D0/                                            GCL1092
C
      ITMAX = NSMAX*100
      IERR = 0
      XABOLD = XAB
      YABOLD = YAB
C   compute XAB, YAB
      CALL LALF0 (BETA,X0,Y0)
      XAB = (1.D0 - ALF) * X0 + ALF * (XSUM + BETA * XDIF)              GCL1092
      YAB = (1.D0 - ALF) * Y0 + ALF * (YSUM + BETA * YDIF)              GCL1092
C
C Set up to do root search for ADIF(S)=0.  ADIF is defined as the cosine
C    of the angle between the normalized gradient vector at s and the 
C    vector from xmep at s to the point on the tunneling path.  For 
C    any point on the concave side of the MEP, ADIF goes to 1 for large
C    negative s, and ADIF goes to -1 for large position s.  Multiple
C    roots of ADIF(s)=0 can exist and we want the one that is the
C    closest to the root for the previous point on the tunneling path.
C
C Evaluate UALF and angles at old s value
      CALL LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
      LCONV = ABS(ADIF) .LT. EPSA .OR. UALF.LT. EPSU
      SSMN = SS(1)
      SSMX = SS(NSMAX)
      DELS = 0.0D0                                                      GCL1092
      IC = 0
C      IDBG = 98
C      WRITE (IDBG,9800) E*CKCAL,ALF,BETA,S,X0,Y0,XAB,YAB
C9800  FORMAT (1X, 131(1H-),/, ' IN LBETAS, E,ALF,BETA=',1P3E13.5,
C     * /, ' S,X0,Y0,XAB,YAB=',/, 1X,1P5E13.5,/,
C     * T15, 'S', T28, 'DELS', T41, 'SSMN', T54, 'SSMX', T67, 'ADIF',
C     * T79, 'UALF', T93, 'LCONV')
C      WRITE (IDBG,9801) IC,S,DELS,SSMN,SSMX,ADIF,UALF,LCONV,X,Y,UX,UY
C9801  FORMAT(1X, I9, 1P6E13.5, 6X, L1, 4E10.2)
      IF (.NOT.LCONV) THEN
C Compute first guess at new s and evaluate UALF and angles there.
C    First guess obtained from the dot product of the vector (DELX,DELY)
C    from the previous and current points on the tunneling path and the 
C    gradient vector at the previous value of s on the MEP.
         DELX = XAB - XABOLD
         DELY = YAB - YABOLD
         DELS = -UY*DELX + UX*DELY
         IF (S+DELS .GT. SSMX) DELS=0.5D0*(SSMX-S)                      GCL1092
         IF (S+DELS .LT. SSMN) DELS=0.5D0*(SSMN-S)                      GCL1092
         DELS = SIGN(MAX(MIN(ABS(DELS),DELSMX),DELSMN),DELS)
         SOLD = S
         ADIFO = ADIF
         S = S + DELS
         UOLD = UALF
C Evaluate UALF and angles at new s
         CALL LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
         LCONV = ABS(ADIF) .LT. EPSA .OR. UALF.LT. EPSU
C         WRITE (IDBG,9801) IC,S,DELS,SSMN,SSMX,ADIF,UALF,LCONV,
C     *      X,Y,UX,UY
C      
      END IF
      IF (.NOT.LCONV .AND. UOLD.GT.DELSMN .AND. UALF.GT.DELSMN) THEN
C Choose sign of DELS so that initial steps along s decrease ADIF
         DELS = -SIGN(ABS(DELS),ADIF*(S-SOLD)/(ADIF-ADIFO))
C Step along s until ADIF changes sign
10       CONTINUE
            IC = IC + 1
            SOLD = S
            ADIFO = ADIF
            UOLD = UALF
15          CONTINUE
               S = S + DELS
               IF (S .GT.SS(NSMAX) .OR. S.LT.SS(1)) THEN
                  S = S - DELS
                  DELS = 0.1D0*DELS                                     GCL1092
                  IF (ABS(DELS).LT.EPSD) THEN
                     WRITE (6, 6000) E*CKCAL,ALF,BETA,S
                     IERR = 1
                     RETURN
                  END IF
               GO TO 15
               END IF
C Evaluate UALF and angles at s
            CALL LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
            LCONV = ABS(ADIF) .LT. EPSA .OR. UALF.LT. EPSU
C            WRITE (IDBG,9801) IC,S,DELS,SSMN,SSMX,ADIF,UALF,LCONV,
C     *       X,Y,UX,UY
            IF (.NOT.LCONV .AND. UALF.GT.DELSMN .AND. 
     *       ADIF*ADIFO .GT. 0.0D0) THEN                                GCL1092
C ADIF did not change sign.  If ADIF is increasing then choose
C   direction of step so that ADIF will eventually go through
C   zero and increase step size to quickly step over the region
C   where ADIF is going through a max.
               IF (ABS(ADIF).GT.ABS(ADIFO)) THEN
                  DELS = SIGN(MIN(2.0D0*ABS(DELS),DELSMX),ADIF)         GCL1092
               END IF
               IF (IC.GT.ITMAX) THEN
                  WRITE (6,*) ' In LBETAS, too many iterations to find',
     *               ' bounds on s, IC,DELS=',IC,DELS
                  IERR = 2
                  RETURN
               END IF
               GO TO 10
            END IF
C ADIF went through zero, set min and max values for root search
         SSMN = MIN(S,SOLD)
         SSMX = MAX(S,SOLD)
      END IF
C
C  Newton-Secant root search for ADIF=0
      IF (.NOT.LCONV .AND. UALF.GT.DELSMN) THEN
30       CONTINUE
            IC = IC + 1
C    Compute Newton-Secant step
            DELS =  -ADIF*(S-SOLD)/(ADIF-ADIFO)
            ADIFO = ADIF
            SOLD = S
            UOLD = UALF
C    Check that Newton-Secant step is not too large and in bounds
            DELS = SIGN(MIN(ABS(DELS),DELSMX),DELS)
            ST = S+DELS
            IF (ST.LE.SSMN .OR. ST.GE.SSMX) DELS = 0.5D0*(SSMN+SSMX)-S  GCL1092
            S = S + DELS
C Evaluate UALF and angles at s
            CALL LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
            LCONV = ABS(ADIF) .LT. EPSA .OR. UALF.LT. EPSU .OR.
     *       ABS(DELS) .LT. EPSD
            IF (ADIFO*ADIF*DELS.LT.0.0D0) THEN                          GCL1092
               SSMX = S
            ELSE
               SSMN = S
            END IF
C            WRITE (IDBG,9801) IC,S,DELS,SSMN,SSMX,ADIF,UALF,LCONV,
C     *       X,Y,UX,UY
            IF (IC.GT.ITMAX) THEN
               WRITE (6,*) ' In LBETAS, too many iterations in Newton-',
     *            'Secant search, IC,DELS=', IC,DELS
               IERR = 2
               RETURN
            END IF
         IF (.NOT.LCONV .AND. UALF.GT.DELSMN) GO TO 30
      END IF
      IF (.NOT.LCONV) THEN
         DELS = SIGN(MIN(DELSMN,ABS(DELS)),DELS)
         IF (UALF.GT.UOLD) DELS = -DELS
C   UALF value too small to compute angle, step along s and
C   minimize ualf
40       CONTINUE
            IC = IC + 1
            S = S + DELS
            UOLD = UALF
C Evaluate UALF and angles at s
            CALL LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
C            WRITE (IDBG,9802) IC,S,DELS,SSMN,SSMX,UALF
C9802        FORMAT(1X, I9, 1P4E13.5, 13X, E13.5)
            LCONV = UALF .LT. EPSU .OR. ABS(DELS) .LT. EPSD
            IF (UALF .GE. UOLD) DELS = -0.5D0 * DELS                    GCL1092
            IF (IC.GT.ITMAX) THEN
               WRITE (6,*) ' In LBETAS, too many iterations in search',
     *            ' for min of UALF, IC,DELS=', IC, DELS
               IERR = 2
               RETURN
            END IF
         IF (.NOT.LCONV) GO TO 40
      END IF
      RETURN
6000  FORMAT(2X,T5,'Warning: In LBETAS, s out of range, E,ALF,BETA,s=',
     *   1P,4E13.5)                                                     GCL0992
      END
