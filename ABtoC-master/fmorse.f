!!!*************************************************************
! 文件/File: fmorse.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: fmorse.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      DOUBLE PRECISION FUNCTION SCALIN(A,AL,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C---------------------------------------
C     NORMALIZING THE MORSE-OSCILLATOR EIGENFUNCTIONS
C     ROUTINE CALCULATES THE NATURAL LOGARITHM OF THE NORMALIZATION
C     FACTOR [MEASURED IN ANGSTROM**(-1/2)]
C----------------------------------------
      INTEGER IR,II,N
      REAL*8 X,XR,RR,Y1,XLNFAC,SC,A,Al,FACLOG
      REAL*8 Y2,GAMMLN

      X=AL+FLOAT(N)+1.0
      XR=X
      Y1=0.0
      IF(X.LT.12.0)GOTO 3
      IR=INT(X/3.0E+00)
      RR=FLOAT(IR)
      XR=X-RR

      DO 2 II=1,IR
  2   Y1=Y1+LOG(X-FLOAT(II))

  3   XLNFAC=FACLOG(N)
      Y2=XLNFAC-GAMMLN(XR)
      SC=LOG(A*AL)+Y2-Y1
      SCALIN=SC/2.0D+00
      RETURN
      END




C --------------------------------------------- C
C --------------------------------------------- C
      REAL*8 FUNCTION FACLOG(N)
      REAL*8 QQ
      INTEGER N,II

      QQ=0.0E+00
      IF (N .LE. 1) THEN
          FACLOG=QQ
          RETURN
      ENDIF
      DO 20 II=2,N
20    QQ = QQ + DLOG(DFLOAT(II))
      FACLOG=QQ
      RETURN
      END
C
C
C
      REAL*8 FUNCTION FMORSE(N,RVAL,RK,R0,A,SCF)
      REAL*8 RVAL,RK,R0,A,SCF,R,X,RN,ALPHA,FCTVAL
      REAL*8 PLAG,POWER_EXP
      INTEGER N
C--------------------------------------------------------
C      CALCULATES MORSE OSCILLATOR EIGENFUNCTIONS
C--------------------------------------------------------
      R=RVAL
      X=2.0E+00*RK*DEXP(-A*(R-R0))
      RN=DFLOAT(N)
      ALPHA=2.0E+00*(RK-RN)-1.0E+00
      POWER_EXP = SCF+(ALPHA*DLOG(X)-X)/2.0E+00
      IF ( POWER_EXP.LT.-100 ) THEN 
          FCTVAL = 0.0 
      ELSE 
         FCTVAL=DEXP(POWER_EXP)*PLAG(N,ALPHA,X)
      ENDIF 
      FMORSE=FCTVAL
      RETURN
      END
C
C
C
      FUNCTION PLAG(N,A,X)
C--------------------------------------------------------
C      CALCULATES ASSOCIATED LAGUERRE POLYNOMIALS
C--------------------------------------------------------
      REAL*8 QQ0,QQ1,A,X,XJJ,QQ2,PLAG
      INTEGER JJ,N

      QQ0=1.0D+00
      IF (N .EQ. 0) THEN
           PLAG=QQ0
           RETURN
      ENDIF
C
      QQ1=A+1.0D+00-X
      IF (N .EQ. 1) THEN
           PLAG=QQ1
           RETURN
      ENDIF
C
      DO 10 JJ=2,N
      XJJ=DFLOAT(JJ)
      QQ2=((2.0D+00*XJJ+A-1.0D+00-X)*QQ1
     1    -(        XJJ+A-1.0D+00  )*QQ0)/XJJ
      QQ0=QQ1
10    QQ1=QQ2
      PLAG=QQ2
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION DLAG(N,A,X)
C--------------------------------------------------------
C      CALCULATES DERIVATIVES (WITH RESPECT TO X)
C      OF ASSOCIATED LAGUERRE POLYNOMIALS
C--------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  QQ0,QQ1,XJJ,QQ2,X,A,PLAG
      INTEGER N,JJ

      QQ0=0.0D+00
      IF (N .EQ. 0) THEN
           DLAG=QQ0
           RETURN
      ENDIF
C
      QQ1=-1.0D+00
      IF (N .EQ. 1) THEN
           DLAG=QQ1
           RETURN
      ENDIF
C
      DO 10 JJ=2,N
      XJJ=DFLOAT(JJ)
      QQ2=((2.0D+00*XJJ+A-1.0D+00-X)*QQ1
     1    -(        XJJ+A-1.0D+00  )*QQ0
     2    -PLAG(JJ-1,A,X))/XJJ
      QQ0=QQ1
10    QQ1=QQ2
      DLAG=QQ2
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION GAMMLN(XX)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER J
      REAL*8   COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
