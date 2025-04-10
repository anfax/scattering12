!!!*************************************************************
! 文件/File: lalfmn.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lalfmn.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LALFMN (IC3D, E, ALFMN, THETMN, THET0, THETH,
     *   THET1, PERC, ILAGTH, IERR, LNONAD,NFINAL)
C
C     LALFMN - find minimum theta(alf)
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C     LALFOR - sort ALF, THET into ALPH, THT arrays.  ALPH is in 
C              increasing order.  Alpha at the minimum theta and the 
C              two closest alphas are kept along with their associated 
C              theta values.  
C     QUADFT - quadratic fit of three points
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LCONV,LNONAD
      DIMENSION ALFMN(2), THETMN(2), THET0(2), THETH(2), THET1(2),
     *PERC(2,2)
      DIMENSION THET(2), PBET(2)
      DIMENSION ALPH(3), THT(3)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /LAGCM2/ TOLLAG
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      SAVE EPSMAX, EPSALF                                               TCA1097
      DATA EPSMAX/1.D-6/, EPSALF/1.D-6/
C
C      IDBG = 96
C      WRITE (IDBG, 900) E*627.5095D0                                   GCL1092
C900   FORMAT (1X,130('-'),/,' LALFMN, E=',1PE15.7,/,4X,'IC',2X,'MIN',
C     *   10X,'ALPH',39X,,'THT')
      DALF = .5D0                                                       GCL1092
      ALF = -.5D0                                                       GCL1092
      IC = 0
      IALFMN = 0
C  loop over alf = 0., .5, and 1., calculate theta and store alf and
C    theta on grid.  IALFMN = index of the minimum theta on the grid
      DO 20 IALF = 1,3
         ALF = ALF + DALF
         CALL LAGTH (E,ALF,THET,IC3D,PBET,IERR,LNONAD,NFINAL)
         ILAGTH = ILAGTH + 1
         THT(IALF) = THET(2)
         ALPH(IALF) = ALF
         IF (IALF .EQ. 1) THEN
            THET0(1) = THET(1)
            THET0(2) = THET(2)
         ELSE IF (IALF .EQ. 2) THEN
            THETH(1) = THET(1)
            THETH(2) = THET(2)
            PERC(1,1) = PBET(1)
            PERC(2,1) = PBET(2)
         ELSE IF (IALF .EQ. 3) THEN
            THET1(1) = THET(1)
            THET1(2) = THET(2)
            PERC(1,2) = PBET(1)
            PERC(2,2) = PBET(2)
         END IF
         IF (IALFMN.EQ.0 .OR. THET(2).LE.THT(IALFMN)) THEN
            IALFMN=IALF
            THETMN(1) = THET(1)
         END IF
20    CONTINUE
C      WRITE (IDBG,910) IC,IALFMN,ALPH,THT
C910   FORMAT (1X, 2I5, 1P6E13.5)
C
      IF (.NOT.LGS2(3)) THEN
C  Try minimum from previous energy
         ALF = ALFMN(2)
         IF (ALF.GT.EPSALF .AND. ABS(ALF-0.5D0).GT.EPSALF .AND.         GCL1092
     *    (1.0D0-ALF).GT.EPSALF) THEN                                   GCL1092
            CALL LAGTH (E,ALF,THET,IC3D,PBET,IERR,LNONAD,NFINAL)
            ILAGTH = ILAGTH + 1
            CALL LALFOR (ALF,THET,ALPH,THT,IALFMN,THETMN)
C            WRITE (IDBG,910) IC,IALFMN,ALPH,THT
         END IF
         IC = 0
30       CONTINUE
            IC = IC + 1
            IF (IALFMN.EQ.1 .OR. (IALFMN.EQ.2 .AND. 
     *         ((ALPH(2)-ALPH(1)).GE.(ALPH(3)-ALPH(2))))) THEN
               ALF = 0.5D0*(ALPH(1)+ALPH(2))                            GCL1092
            ELSE
               ALF = 0.5D0*(ALPH(2)+ALPH(3))                            GCL1092
            END IF
            CALL LAGTH (E,ALF,THET,IC3D,PBET,IERR,LNONAD,NFINAL)
            ILAGTH = ILAGTH + 1
C  store alpha at minimum theta and two closest alphas and theta(alpha)
            CALL LALFOR (ALF,THET,ALPH,THT,IALFMN,THETMN)
C            WRITE (IDBG,910) IC,IALFMN,ALPH,THT
            LCONV = (ALPH(3)-ALPH(1)).LE.EPSALF
            IF (.NOT.LCONV) THEN
C  evaluate the convergence criterion, EPS
               T = EXP(-2.0D0*THT(IALFMN))                              GCL1092
               EPS = 0.01D0                                             GCL1092
               IF(T .GT. TOLLAG) EPS=-0.5D0*LOG(1.0D0-TOLLAG/T)/        GCL1092
     1                                THT(IALFMN)                       GCL1092
               EPS = MAX(EPS, EPSMAX)
C   check for convergence in theta: if quadratic fit to theta differs
C   from computed theta less than eps or change in theta from last
C   minimum theta is less than eps
               LCONV = ABS(1.0D0-THT(2)/THT(1)).LE.EPS .AND.            GCL1092
     *                 ABS(1.0D0-THT(3)/THT(2)).LE.EPS                  GCL1092
            END IF
         IF (.NOT.LCONV .AND. IC.LT.50) GO TO 30
         IF (.NOT.LCONV) WRITE (6, 6000) ALPH, THT
      END IF
      ALFMN(1) = ALPH(IALFMN)
      THETMN(2) = THT(IALFMN)
      ALFMN(2) = ALPH(IALFMN)
      RETURN
 6000 FORMAT(2X,T5,'Warning: LAG not converged with respect to alpha',
     1  /' ALF=',1P,3E13.5,' THETA=',3E13.5)                            GCL0992
      END
