!!!*************************************************************
! 文件/File: clebsh.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: clebsh.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      THIS FILE INCLUDES SUBROUTINES DEVOTED TO CALCULATION OF THE C
C      MATRIX ELEMENTS ON THE ROTATION   BASE FUNCTIONS             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATES MATRIX ELEMENTS OF THE irreducible repres. D(1)  
C             : of the 3D rot. group 
C             : ON THE ROTATION FUNCTIONS <J K M tau |D_{ss'}^{1} | J' K' M' tau'>
C             : 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE Calc_RotateMatrElement(ThreeJsymbol,SixJsymbol)

      include 'dipolsys.h'
      include 'molcul.h'

      INTEGER N1,K1,K2,N2,J2,J1

      INTEGER DeltaJ,DeltaN,DeltaK

      REAL*8 ThreeJsymbol(KMAXAB+1,3,KMAXAB+1,3),
     1       SixJsymbol(KMAXAB+1,KMAXAB+1,3,3),
     1       THREEJ,SIXJ
C
      REAL*8 F_kk,ProdN1N2,Spin,J_1,J_2


C =======================================================C
C ================ START SUBROUTINE =====================C




      DO  DeltaJ = -1,1,1
      DO  DeltaN = -1,1,1        
      DO  J2 = 0,KMAXAB 
      DO  N2 = 0,KMAXAB

          Spin = (XMULTI-1.0D+00)*0.5D+00

          IF ( MOD(INT(XMULTI),2).EQ.0) THEN 
               J_1 = DFLOAT(J2+DeltaJ)+0.5D+00
               J_2 = DFLOAT(J2)+0.5D+00
          ELSE 
               J_1 = DFLOAT(J2+DeltaJ)
               J_2 = DFLOAT(J2)
          ENDIF 

          N1 = N2 + DeltaN

          IF ( J2+DeltaJ.GE.0 .AND. N2+DeltaN.GE.0  ) THEN 

            ProdN1N2 = DSQRT( (2.0D+00*DFLOAT(N1)+1.0D+00)*
     1                        (2.0D+00*DFLOAT(N2)+1.0D+00) )

            SixJsymbol(J2+1,N2+1,DeltaN+2,DeltaJ+2) = ProdN1N2 *
     1                            SIXJ(J_2, DFLOAT(N2) , Spin ,
     1                                 DFLOAT(N1) , J_1, 1.0D+00)

          ELSE
            SixJsymbol(J2+1,N2+1,DeltaN+2,DeltaJ+2) = 0.0
          ENDIF 


      ENDDO     ! --- N1
      ENDDO     ! --- J1
      ENDDO     ! --- N2
      ENDDO     ! --- J2



      DO  DeltaK = -1,1,1
      DO  DeltaN = -1,1,1
      DO  N2 = 0,KMAXAB 
      DO  K2 = 0,KMAXAB

          N1 = N2 + DeltaN
          K1 = K2 + DeltaK

          IF ( N1.GE.0 .AND. K1.GE.0  ) THEN

          F_kk = 1.0E+00
          IF (K2+K1-1.EQ.0) F_kk = dsqrt(2.0D+00)

          ThreeJsymbol(N2+1,DeltaN+2,K2+1,DeltaK+2) = F_kk*
     1                        THREEJ(N2,1,N1,K2,DeltaK,-K1) 

          ENDIF
     
      ENDDO     ! --- N1
      ENDDO     ! --- J1
      ENDDO     ! --- N2
      ENDDO     ! --- J2


C
      RETURN
      END
C###################        RotateMatrElement      #########################
C###########################################################################



C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : 
C UPDATES     :
C LANGUAGE    : FORTRAN
C             : CALCULATES WIGNER 3-J SYMBOLS USING FORMULA OF P34 OF BRINK &
C             : SATCHLER IN "ANGULAR MOMENTUM"
C             : NOTE FAKT(A)=FACT(A)/10**A
C             : A=J1,B=J2,C=J,AL=M1,BE=M2,GA=-M
C             : 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      REAL*8 FUNCTION THREEJ(J1,J2,J3,K1,K2,K3)
C      DOUBLE PRECISION FUNCTION THREEJ(A,B,C,AL,BE,GA)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER J1,J2,J3,K1,K2,K3,NEWMIN,NEWMAX,NEW,IPHASE
      REAL*8 A,B,C,AL,BE,GA,DELTA,FAKT,CLEBSH,Minus
      REAL*8 TERM,TERM1,TERM2,TERM3,SUMM,DNEW,TERM4,TERM5,TERM6

C
      A = J1
      B = J2
      C = J3
      AL= K1
      BE= K2
      GA= K3

      THREEJ=0.D0
C
C     (J1+J2).GE.J AND J.GE.ABS(A-B)    -M=M1+M2    J1,J2,J.GE.0
C     ABS(M1).LE.J1    ABS(M2).LE.J2   ABS(M).LE.J
C
      IF(C.GT.A+B)RETURN
      IF(C.LT.DABS(A-B))RETURN
      IF(A.LT.0.D0.OR.B.LT.0.D0.OR.C.LT.0.D0)RETURN
      IF(A.LT.DABS(AL).OR.B.LT.DABS(BE).OR.C.LT.DABS(GA))RETURN
      IF(-1.D0*GA.NE.AL+BE)RETURN
C
C
C     COMPUTE DELTA(ABC)
C
      DELTA=DSQRT(FAKT(A+B-C)*FAKT(A+C-B)*FAKT(B+C-A)/FAKT(A+B+C+1.D0))
C
C
      TERM1=FAKT(A+AL)*FAKT(A-AL)
      TERM2=FAKT(B-BE)*FAKT(B+BE)
      TERM3=FAKT(C+GA)*FAKT(C-GA)
      TERM=DSQRT((2.D0*C+1.D0)*TERM1*TERM2*TERM3)
C
C
C     NOW COMPUTE SUMMATION TERM
C
C     SUM TO GET SUMMATION IN EQ(2.34) OF BRINK AND SATCHLER.  SUM UNTIL
C     A TERM INSIDE FACTORIAL GOES NEGATIVE.  NEW IS INDEX FOR SUMMATION
C     .  NOW FIND WHAT THE RANGE OF NEW IS.
C
C
      NEWMIN=IDNINT(DMAX1((A+BE-C),(B-C-AL),0.D0))
      NEWMAX=IDNINT(DMIN1((A-AL),(B+BE),(A+B-C)))
C
C
      SUMM=0.D0+00
C
C
             DO 10 NEW=NEWMIN,NEWMAX
             DNEW=DFLOAT(NEW)
             TERM4=FAKT(A-AL-DNEW)*FAKT(C-B+AL+DNEW)
             TERM5=FAKT(B+BE-DNEW)*FAKT(C-A-BE+DNEW)
             TERM6=FAKT(DNEW)*FAKT(A+B-C-DNEW)
             SUMM=SUMM+(-1.D0)**NEW/(TERM4*TERM5*TERM6)
10           CONTINUE
C
C     SO CLEBSCH-GORDON <J1J2M1M2LJM> IS CLEBSH
C
      CLEBSH=DELTA*TERM*SUMM/DSQRT(10.D0)
C
C     CONVERT CLEBSCH-GORDON TO THREEJ
C
      IPHASE=IDNINT(A-B-GA)
      Minus = -1.0D+00 
      IF (Mod(IPHASE,2).EQ.0) Minus = 1.0D+00
      THREEJ=Minus*CLEBSH/DSQRT(2.D0*C+1.D0)


C     THREEJ=(-1.D0)**(IPHASE)*CLEBSH/DSQRT(2.D0*C+1.D0)
C
C
      RETURN
      END
C######################          THREEJ           ##########################
C###########################################################################

C
C

C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C             : THIS FUNCTION COMPUTES N!/(10**N).  STOP 101 IMPLIES TRIED
C             : TO TAKE FACTORIAL OF A NEGATIVE NO.
C             : 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################


      REAL*8 FUNCTION FAKT(A)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A,AX
      INTEGER I,IC
C

C
      AX=A
      FAKT=1.D0
      IF(AX.EQ.0.D0)RETURN
      FAKT=.1D0
      IF(AX.LT.0.D0) GOTO 9000
C
      IC=IDNINT(AX)
      AX=AX/10.D0
      FAKT=AX
             DO 10  I=1,IC-1
             FAKT=FAKT*(AX-DFLOAT(I)*.1D0)
10           CONTINUE
      RETURN
9000  WRITE (*,9014) AX
9014  FORMAT(1H0,' RENNER.FKT.ERR  NEGATIVE ARGUMENT FOR FUNCTI',
     1      'ON FAKT. ARGUMENT = ',E12.5)
      STOP
      END


C######################      FAKT                 ##########################
C###########################################################################




C##########################################################################
C             :
C LANGUAGE    : FORTRAN
C     ROUTINE CALCULATES  6J SYMBOLS {X1 X2 X3}
C                                    {X4 X5 X6}
C
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      DOUBLE PRECISION FUNCTION SIXJ(X1,X2,X3,X4,X5,X6)
      REAL*8 AA1,AA2,AA3,AA4,X1,X2,X3,X4,X5,X6,FRAC,FACLOG,TERM,
     1       EXPON,SUM,XK,YK
      INTEGER I1,I2,I3,I4,K1,K2,K3,KMIN,KMAX,
     1        I2X1,I2X2,I2X3,I2X4,I2X5,I2X6,K
C
      SIXJ=0.D0

C
      AA1=(-X1+X2+X3)*(X1-X2+X3)*(X1+X2-X3)
      IF (AA1.LT.0.D0) GO TO 1000
      AA2=(-X1+X5+X6)*(X1-X5+X6)*(X1+X5-X6)
      IF (AA2.LT.0.D0) GO TO 1000
      AA3=(-X4+X2+X6)*(X4-X2+X6)*(X4+X2-X6)
      IF (AA3.LT.0.D0) GO TO 1000
      AA4=(-X4+X5+X3)*(X4-X5+X3)*(X4+X5-X3)
      IF (AA4.LT.0.D0) GO TO 1000
C
C
      I1=IDINT((X1+X2+X3))
      FRAC=X1+X2+X3-DFLOAT(I1)
      IF (DABS(FRAC).GT.1.D-9) GO TO 1000
      I2=IDINT((X1+X5+X6))
      FRAC=X1+X5+X6-DFLOAT(I2)
      IF (DABS(FRAC).GT.1.D-9) GO TO 1000
      I3=IDINT((X4+X2+X6))
      FRAC=X4+X2+X6-DFLOAT(I3)
      IF (DABS(FRAC).GT.1.D-9) GO TO 1000
      I4=IDINT((X4+X5+X3))
      FRAC=X4+X5+X3-DFLOAT(I4)
      IF (DABS(FRAC).GT.1.D-9) GO TO 1000
C
      K1=IDINT((X2+X3+X5+X6))
      K2=IDINT((X1+X3+X4+X6))
      K3=IDINT((X1+X2+X4+X5))
C
      KMIN=MAX0(I1,I2,I3,I4)
      KMAX=MIN0(K1,K2,K3)
C
      IF (KMIN.GT.KMAX) GO TO 1000
C
      I2X1=IDINT((2.D0*X1+0.1D0))
      I2X2=IDINT((2.D0*X2+0.1D0))
      I2X3=IDINT((2.D0*X3+0.1D0))
      I2X4=IDINT((2.D0*X4+0.1D0))
      I2X5=IDINT((2.D0*X5+0.1D0))
      I2X6=IDINT((2.D0*X6+0.1D0))
C
C
      EXPON=0.5D0*(FACLOG(I1-I2X1)+FACLOG(I1-I2X2)+FACLOG(I1-I2X3)
     +    +FACLOG(I2-I2X1)+FACLOG(I2-I2X5)+FACLOG(I2-I2X6)
     +    +FACLOG(I3-I2X4)+FACLOG(I3-I2X2)+FACLOG(I3-I2X6)
     +    +FACLOG(I4-I2X4)+FACLOG(I4-I2X5)+FACLOG(I4-I2X3)
     +    -FACLOG(I1+1)-FACLOG(I2+1)-FACLOG(I3+1)-FACLOG(I4+1))
     +    +FACLOG(KMIN+1)
     +    -FACLOG(KMIN-I1)-FACLOG(KMIN-I2)-FACLOG(KMIN-I3)
     +    -FACLOG(KMIN-I4)
     +    -FACLOG(K1-KMIN)-FACLOG(K2-KMIN)-FACLOG(K3-KMIN)
      TERM=DFLOAT((-1)**KMIN)*DEXP(EXPON)
C
C
      SIXJ=TERM
C
      IF (KMAX.EQ.KMIN) GO TO 1000
C
      SUM=1.D0
      TERM=1.D0
      DO 200 K=(KMIN+1),KMAX
      XK=DFLOAT((K+1)*(K1-K+1)*(K2-K+1)*(K3-K+1))
      YK=DFLOAT((K-I1)*(K-I2)*(K-I3)*(K-I4))
      TERM=-TERM*XK/YK
      SUM=SUM+TERM
  200 CONTINUE
C
      SIXJ=SIXJ*SUM
1000  RETURN
C
  900 WRITE(6,6111)
 6111 FORMAT(1H0,' RENNER.SXJ.ERR   SUM OF J VALUES IN A TRIANGLE ',
     1 'OF 6J NON-INTEGRAL ')
      WRITE(6,6112) X1,X2,X3,X4,X5,X6
 6112 FORMAT(1H0,'                  J VALUES:', 6F10.3)
C
c      STOP
      END



C
      DOUBLE PRECISION FUNCTION THREEJFLOAT(A,B,C,AL,BE,GA)
C
C     CALCULATES WIGNER 3-J SYMBOLS USING FORMULA OF P34 OF BRINK &
C     SATCHLER IN "ANGULAR MOMENTUM"
C     NOTE FAKT(A)=FACT(A)/10**A
C     A=J1,B=J2,C=J,AL=M1,BE=M2,GA=-M
C
C
      REAL*8 A,B,C,AL,BE,GA,DELTA,TERM1,TERM2,TERM3,TERM,
     1       CLEBSH,IPHASE,FAKT,NEWMIN,NEWMAX,SUMM,DNEW,
     1       TERM4,TERM5,TERM6
      INTEGER NEW

C
      THREEJFLOAT=0.D0
C
C     (J1+J2).GE.J AND J.GE.ABS(A-B)    -M=M1+M2    J1,J2,J.GE.0
C     ABS(M1).LE.J1    ABS(M2).LE.J2   ABS(M).LE.J
C
      IF(C.GT.A+B)RETURN
      IF(C.LT.DABS(A-B))RETURN
      IF(A.LT.0.D0.OR.B.LT.0.D0.OR.C.LT.0.D0)RETURN
      IF(A.LT.DABS(AL).OR.B.LT.DABS(BE).OR.C.LT.DABS(GA))RETURN
      IF(-1.D0*GA.NE.AL+BE)RETURN
C
C
C     COMPUTE DELTA(ABC)
C
      DELTA=DSQRT(FAKT(A+B-C)*FAKT(A+C-B)*FAKT(B+C-A)/FAKT(A+B+C+1.D0))
C
C
      TERM1=FAKT(A+AL)*FAKT(A-AL)
      TERM2=FAKT(B-BE)*FAKT(B+BE)
      TERM3=FAKT(C+GA)*FAKT(C-GA)
      TERM=DSQRT((2.D0*C+1.D0)*TERM1*TERM2*TERM3)
C
C
C     NOW COMPUTE SUMMATION TERM
C
C     SUM TO GET SUMMATION IN EQ(2.34) OF BRINK AND SATCHLER.  SUM UNTIL
C     A TERM INSIDE FACTORIAL GOES NEGATIVE.  NEW IS INDEX FOR SUMMATION
C     .  NOW FIND WHAT THE RANGE OF NEW IS.
C
C
      NEWMIN=IDNINT(DMAX1((A+BE-C),(B-C-AL),0.D0))
      NEWMAX=IDNINT(DMIN1((A-AL),(B+BE),(A+B-C)))
C
C
      SUMM=0.D0
C
C
             DO 10 NEW=NEWMIN,NEWMAX
             DNEW=DFLOAT(NEW)
             TERM4=FAKT(A-AL-DNEW)*FAKT(C-B+AL+DNEW)
             TERM5=FAKT(B+BE-DNEW)*FAKT(C-A-BE+DNEW)
             TERM6=FAKT(DNEW)*FAKT(A+B-C-DNEW)
             SUMM=SUMM+(-1.D0)**NEW/(TERM4*TERM5*TERM6)
10           CONTINUE
C
C     SO CLEBSCH-GORDON <J1J2M1M2LJM> IS CLEBSH
C
      CLEBSH=DELTA*TERM*SUMM/DSQRT(10.D0)
C
C     CONVERT CLEBSCH-GORDON TO THREEJ
C
      IPHASE=IDNINT(A-B-GA)
      THREEJFLOAT=(-1.D0)**(IPHASE)*CLEBSH/DSQRT(2.D0*C+1.D0)
C
C
      RETURN
      END
