!!!*************************************************************
! 文件/File: dbwlrt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: dbwlrt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE DBWLRT (E2, XN, A, B, XMU, ITYP, EMIN, EMAX, DEL, SGN)
C
C     DBWLRT - performs root search for semiclassical eigenvalue of
C        harm.-quart.  potential.
C
C  Called by:
C     DBWL   - compute semiclassical eigenvalue of harm.-quart.
C        potential.
C
C  Calls:
C     PHID   - compute phi function needed for uniform semiclassical
C     THETA2 - phase integral for semiclassical eigenvalues of double
C              well
C     THETA3 - phase integral for semiclassical eigenvalues of double
C              well
C
C  ITYP = 2 for FB>0
C  ITYP = 3 for FB<0
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSET
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      SAVE EPS                                                          TCA1097
      DATA EPS/1.D-8/
C
      IC = 0
   10 CONTINUE
         IC = IC + 1
         IF (IC .GT. 50) THEN
            WRITE (6, 601) IC, E2, FN2
            RETURN
         END IF
         IF (ITYP .EQ. 2) CALL THETA2(XMU, A, B, E2, TH1, DTH1, TH2,
     *      DTH2)
         IF (ITYP .EQ. 3) CALL THETA3(XMU, A, B, E2, TH1, DTH1, TH2,
     *      DTH2)
         CALL PHID(TH2, PH, DPH)
         IF (TH2 .GT. 27.0D0) THEN
            EX = 1.D35
            ATNEX = 0.5D0*PI
         ELSE
            EX = EXP(PI*TH2)
            ATNEX = ATAN(EX)
         END IF
         FN2 = TH1 - XN - (PH + SGN*ATNEX)/TPI
         IF (ABS(FN2) .LT. EPS) RETURN
         IF (EX .GT. 1.D18) THEN
            DATNEX = 0.5D0/EX
         ELSE
            DATNEX = 0.5D0*EX/(1.0D0+EX*EX)
         END IF
         DFN = DTH1 - (DPH/TPI + SGN*DATNEX)*DTH2
         LSET = .FALSE.
         IF (IC .GT. 1) THEN
C  fit F(E) to quadratic in E using F1,F2,DFN2
            T = E1 - E2
            IF (ABS(T) .GT. 1.D-12) THEN
               C = ((FN1-FN2)/T-DFN)/T
               IF (ABS(C) .GT. 1.D-8) THEN
                  T1 = 0.5D0*DFN/C
                  T = 1.0D0-FN2/(T1*T1*C)
                  IF (T .GE. 0.0D0) THEN
                     E1 = E2
                     FN1 = FN2
                     E2 = E2 - T1*(1.0D0-SQRT(T))
                     LSET = .TRUE.
                  END IF
               END IF
            END IF
         END IF
         IF (.NOT.LSET) THEN
            E1 = E2
            FN1 = FN2
            E2 = E2 - FN2/DFN
         END IF
         IF (E2 .LE. 0.D0) THEN
            DEL = 0.1D0*DEL
            IF (DEL .LT. 1.D-12) THEN
               WRITE (6, 600)
               E2 = 0.0D0
               RETURN
            END IF
            E2 = DEL
         END IF
         IF (FN1 .GT. 0.D0) EMAX = MIN(E1, EMAX)
         IF (FN1 .LT. 0.D0) EMIN = MAX(E1, EMIN)
         IF (E2 .LE. EMIN .OR. E2 .GE. EMAX)
     *      E2 = 0.5D0*(EMIN+EMAX)
      GO TO 10
C
600   FORMAT(/,2X,T5,'Warning: In DBWLRT there is a problem with the ',
     *               'root search for the eigenvalue',
     *       /,2X,T14,'The eigenvalue has been set equal to zero')
601   FORMAT(/,2X,T5,'Warning: In DBWLRT after ',I3,
     *               'iterations: E=',1PE13.5, ' F(E)=', 1PE13.5)
      END
