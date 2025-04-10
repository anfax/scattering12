!!!*************************************************************
! 文件/File: curvex.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: curvex.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE CURVEX
C
C     CURVEX - exponential extrapolation of curvature to asymptotic
C        regions
C
C  Called by:
C     RPHSUM - summarize reaction path info
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      COMMON /CURVE2/ SM2, SP2
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL1096
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
C  interpolate curvature at saddle point
      IF (NSAD .NE. 0) THEN
         DO 10 I = 1,NSAD
            IS = NSSP(I)
            CAPS(IS) = CAPS(IS-1)
            IF (LSYM .AND. IS .EQ. NMID) GO TO 10
            CAPS(IS) = CAPS(IS-1) + (CAPS(IS+1)-CAPS(IS-1))*
     *         (SS(IS)-SS(IS-1))/(SS(IS+1)-SS(IS-1))
   10    CONTINUE
      END IF
C  exponential extrapolation on reactant side
      IF (SM .GT. SM2) THEN
         NSHLF = NSMAX
         IF (NSAD .GE. 1)  NSHLF = NSSP(1) - 1
         DO 20 IS = 1,NSHLF
            ISL = IS
            IF (CAPS(IS) .NE. 0.D0) GO TO 30
   20    CONTINUE
   30    CONTINUE
         IF (CAPS(ISL) .NE. 0.D0 .AND. ISL .GT. 1) THEN                 GCL1092
            ALF = LOG(2.0D0)/(SM-SM2)
            XK0 = CAPS(ISL)*EXP(-ALF*SS(ISL))
            ISL = ISL - 1
            WRITE(6,600) SS(1),SS(ISL),SM2,SM
            DO 40 IS = 1,ISL
               CAPS(IS) = XK0*EXP(ALF*SS(IS))
   40       CONTINUE
         END IF
      END IF
      IF (.NOT. LSYM .AND. SP .LT. SP2) THEN
C  exponential extrapolation on product side
         NSHLF = 1
         IF (NSAD .GT. 0) NSHLF = NSSP(NSAD) + 1
         ISL = NSMAX + 1
         DO 50 IS = NSHLF,NSMAX
            ISL = ISL - 1
            IF (CAPS(ISL) .NE. 0.D0) GO TO 60
   50    CONTINUE
   60    CONTINUE
         IF (CAPS(ISL) .NE. 0.D0 .AND. ISL .LT. NSMAX) THEN             GCL1092
            ALF = LOG(2.0D0)/(SP-SP2)
            XK0 = CAPS(ISL)*EXP(-ALF*SS(ISL))
            ISL = ISL + 1
            WRITE(6,600) SS(ISL),SS(NSMAX),SP2,SP
            DO 70 IS = ISL,NSMAX
               CAPS(IS) = XK0*EXP(ALF*SS(IS))
   70       CONTINUE
         END IF
      END IF
      RETURN
600   FORMAT(/,2X,T5,'The curvature is extrapolated from s = ',F9.4,
     *               'to s = ',F9.4,
     *       /,2X,T5,'by fitting to an exponential, such that the ',
     *       'curvature at s = ',F9.4,' is half that at s = ',F9.4) 
      END
