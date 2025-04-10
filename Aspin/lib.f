!!!*************************************************************
! 文件/File: lib.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lib.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE INVLEG(N,PLEG1,GAM)
C     CALCULATES GAUSS-LEGENDRE QUADRATURE WEIGHTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/HOMO/IHOMO,ISYMM
      DIMENSION GAM(N)
      DIMENSION P(2*N)
      DIMENSION IPIV(N),WORK(3*N)
      DIMENSION PLEG1(N,N)
      LOGICAL IHOMO
      LWORK=3*N
      NP=N
      IF(IHOMO) NP=2*N-1
      DO I=1,N
         TET=GAM(I)
         CALL PLEG0(TET,NP,P)
         DO  J=1,N
            JJ=J-1
            IF(IHOMO)JJ=2*JJ
            PLEG1(I,J)=P(JJ+1)
         ENDDO
      ENDDO
      CALL DGETRF(N,N,PLEG1,N,IPIV,INFO)
      CALL DGETRI(N,PLEG1,N,IPIV,WORK,LWORK,INFO)
      IF (INFO.GT.0) STOP 'INV FAILED'
      RETURN
      END     

C............................................................
 
      SUBROUTINE PLEG0(TET,N,P)
C     ORTHOGONALIZED LEGENDRE POLYNOMIALS
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(N+1)
      X=DCOS(TET)
      P(1)=1.D0
      P(2)=X
      DO KS=3,N
         K=KS-1
         A1=DFLOAT(2*K-1)
         A2=DFLOAT(K-1)
         A3=DFLOAT(K)
         P(KS)=(A1*P(K)*X-A2*P(K-1))/A3
      ENDDO
      RETURN
      END

C     .................................................................
 
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),Y2(N),U(N)
      IF (YP1.GT..99D30) THEN
         Y2(1)=0.D0
         U(1)=0.D0
      ELSE
         Y2(1)=-0.5D0
         U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.D0
         Y2(I)=(SIG-1.D0)/P
         U(I)=(6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
 11   CONTINUE
      IF (YPN.GT..99D30) THEN
         QN=0.D0
         UN=0.D0
      ELSE
         QN=0.5D0
         UN=(3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
 12   CONTINUE
      RETURN
      END

C     .................................................................

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
 1    IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
         GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.D0) PAUSE 'BAD XA INPUT.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.D0
      RETURN
      END

C     .................................................................

      SUBROUTINE DSYFIL(UPLO, N, A, LDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 UPLO
      DIMENSION A(LDA,N)
C
C     SUBROUTINE TO FILL IN THE SECOND TRIANGLE OF A SYMMETRIC MATRIX.
C     IF UPLO='L', THE LOWER TRIANGLE IS FILLED IN
C     IF UPLO='U', THE UPPER TRIANGLE IS FILLED IN
C
      IF(UPLO.EQ.'L') THEN
         DO J=1,N-1
            CALL DCOPY(N-J,A(J,J+1),LDA,A(J+1,J),1)
         ENDDO
      ELSEIF(UPLO.EQ.'U') THEN
         DO J=1,N-1
            CALL DCOPY(N-J,A(J+1,J),1,A(J,J+1),LDA)
         ENDDO
      ENDIF
C
      RETURN
      END

C     .................................................................

      subroutine syminv(a,ia,n)
C     SIMULATES SYMINV SYMMETRIC MATRIX INVERTER WITH LAPACK CALLS
C     THIS VERSION USES ONLY THE UPPER TRIANGLE OF A
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(LWORK=100000000)
      DIMENSION A(IA,N)
      DIMENSION WORK(LWORK),IPIV(N)
      NB=ILAENV(1,'DSYTRF','L',N,-1,-1,-1)
      LWREQ=N*NB
      IF(LWORK.LT.LWREQ) THEN
         WRITE(6,100) LWORK,N,NB
 100     FORMAT(' *** ERROR: ONLY',I8,' WORDS OF WORKSPACE AVAILABLE',
     1       ' IN SYMINV.'/'LAPACK ROUTINE DSYTRF NEEDS AT LEAST N*NB:',
     2       ' N =',I5,' AND NB =',I5,' ON THIS CALL.')
         STOP
      ENDIF
      CALL DSYTRF('L',N,A,IA,IPIV,WORK,LWORK,INFO)
      IF (INFO .NE. 0) THEN
         WRITE (6,120) INFO
 120     FORMAT(' *** ERROR IN DSYTRF: INFO =',I3)
         STOP
      END IF
      CALL DSYTRI('L',N,A,IA,IPIV,WORK,INFO)
      IF (INFO .NE. 0) THEN
         WRITE (6,130) INFO
 130     FORMAT(' *** ERROR IN DSYTRI: INFO =',I3)
         STOP
      END IF
      RETURN
      END

C............................................................

      SUBROUTINE KSYM(AK,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK(N,N)
      DO I=1,N
         DO J=1,I
            TMP=0.5D0*(AK(I,J)+AK(J,I))
            AK(I,J)=TMP
            AK(J,I)=TMP
         ENDDO
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------

      SUBROUTINE FACTORIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/FCT/FACT(100000)
      FACT(1)=0.D0
      FACT(2)=0.D0
      N=10000
 2    DO I=3,N
         FACT(I)=FACT(I-1)+DLOG(DFLOAT(I-1))
      ENDDO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION XF6J(A,B,E,D,C,F)
*
*                                    | A  B  E |
*   PROGRAM TO COMPUTE THE 6J SYMBOL {         }
*                                    | D  C  F |
*   AUTHOR: B. FOLLMEG
*   CURRENT REVISION DATE: 4-MAY-1997
*
* -------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FCT/ SI(100000)
      DATA TOL,ZERO,ONE,TWO /1.D-10,0.D0,1.D0,2.D0/
      X=ZERO
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( A B E)
      IF ((E .GT. (A + B)) .OR. (E .LT. ABS(A - B))) GOTO 40
      SUM = A + B + E
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IABE = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( D C E)
      IF ((E .GT. (C + D)) .OR. (E .LT. ABS(C - D))) GOTO 40
      SUM = D + C + E
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IDCE = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( A C F)
      IF ((F .GT. (A + C)) .OR. (F .LT. ABS(A - C))) GOTO 40
      SUM = A + C + F
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IACF = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( D B F)
      IF ((F .GT. (D + B)) .OR. (F .LT. ABS(D - B))) GOTO 40
      SUM = D + B + F
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IDBF = NINT(SUM)
      IABDC = NINT(A + B + D + C)
      IAEDF = NINT(A + E + D + F)
      IBECF = NINT(B + E + C + F)
      MINCHI = MAX(IABE,IDCE,IACF,IDBF,-1)
      MAXCHI = MIN(IABDC,IAEDF,IBECF) - MINCHI
* INDICES FOR DELTAS
      DELTA = ZERO
      I2A = NINT(TWO*A) - 1
      I2B = NINT(TWO*B) - 1
      I2C = NINT(TWO*C) - 1
      I2D = NINT(TWO*D) - 1
      I2E = NINT(TWO*E) - 1
      I2F = NINT(TWO*F) - 1
* DELTA(ABE)
      J1 = IABE - I2A
      J2 = IABE - I2B
      J3 = IABE - I2E
      J4 = IABE + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(DCE)
      J1 = IDCE - I2D
      J2 = IDCE - I2C
      J3 = IDCE - I2E
      J4 = IDCE + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(ACF)
      J1 = IACF - I2A
      J2 = IACF - I2C
      J3 = IACF - I2F
      J4 = IACF + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(DBF)
      J1 = IDBF - I2D
      J2 = IDBF - I2B
      J3 = IDBF - I2F
      J4 = IDBF + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
      DELTA = 0.5D0 * DELTA
      IABDC = IABDC - MINCHI
      IAEDF = IAEDF - MINCHI
      IBECF = IBECF - MINCHI
      IABE = MINCHI - IABE
      IDCE = MINCHI - IDCE
      IACF = MINCHI - IACF
      IDBF = MINCHI - IDBF
      ABDC = DBLE(IABDC - MAXCHI)
      AEDF = DBLE(IAEDF - MAXCHI)
      BECF = DBLE(IBECF - MAXCHI)
      ABE = DBLE(MAXCHI + IABE + 1)
      DCE = DBLE(MAXCHI + IDCE + 1)
      ACF = DBLE(MAXCHI + IACF + 1)
      DBF = DBLE(MAXCHI + IDBF + 1)
* LOOP OVER CHI
      X = 1.D0
      IPOWER = 0
      IF (MAXCHI .LE. 0) GOTO 30
      II = MINCHI + MAXCHI + 2
      DO 10 ICHI = 1, MAXCHI
         XCHI = DBLE(ICHI)
         XA = (ABDC + XCHI) * (AEDF + XCHI) * (BECF + XCHI)
         XB = (ABE - XCHI) * (DCE - XCHI) * (ACF - XCHI) * (DBF - XCHI)
         X = 1.D0 - XA * (II - ICHI) * X / XB
10    CONTINUE
      IF (X) 20,40,30
20    X = -X
      IPOWER = 1
30    X = DLOG(X) + SI(MINCHI+2)-SI(IABDC+1)-SI(IAEDF+1)-SI(IBECF+1)
     :            - SI(IABE+1)-SI(IDCE+1)-SI(IACF+1)-SI(IDBF+1) + DELTA
      IPOWER = IPOWER + MINCHI
      X = (-1.D0)**IPOWER * DEXP(X)
40    XF6J = X
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION XF3J(A,B,C,AM,BM,CM)
*     PROGRAM TO COMPUTE THE 3J SYMBOL (A B C / AM BM CM)
*     AUTHORS:  T. ORLIKOWSKI AND B. FOLLMEG
*     CURRENT REVISION DATE: 4-MAY-1997

* --------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FCT/ SI(100000)
      DATA TOL /1.D-10/
      DATA ZERO,ONE,TWO /0.D0, 1.D0,2.D0/
      X = ZERO
* CHECK FOR TRIANGULAR CONDITIONS
      IF ((C .GT. (A + B)) .OR. (C .LT. ABS(A - B))) GOTO 3
      ABC = A + B + C
      IABC = NINT(ABC)
      IF (ABS(AM + BM + CM) .GT. TOL) GOTO 3
      IF (ABS(AM) .GT. A) GOTO 3
      IF (ABS(BM) .GT. B) GOTO 3
      IF (ABS(CM) .GT. C) GOTO 3
      IF (AM.EQ.0 .AND. BM.EQ.0 .AND.((IABC/2)*2).NE.IABC) GOTO 3
      IACBM = NINT(A - C + BM)
      IBCAM = NINT(B - C - AM)
      IABMC = NINT(A + B - C)
      IAMAM = NINT(A - AM)
      IBPBM = NINT(B + BM)
      IAPAM = NINT(A + AM) + 1
      IBMBM = NINT(B - BM) + 1
      ICPCM = NINT(C + CM) + 1
      ICMCM = NINT(C - CM) + 1
      MINCHI = MAX0(0,IBCAM,IACBM)
      MAXCHI = MIN0(IABMC,IAMAM,IBPBM) - MINCHI
      IABMC = IABMC - MINCHI
      IAAM = IAMAM - MINCHI
      IBBM = IBPBM - MINCHI
      IBCAM = MINCHI - IBCAM
      IACBM = MINCHI - IACBM
      IAMAM = IAMAM + 1
      IBPBM = IBPBM + 1
* COMPUTE DELTA
      J1 = IABC - NINT(TWO*A) + 1
      J2 = IABC - NINT(TWO*B) + 1
      J3 = IABC - NINT(TWO*C) + 1
      J4 = IABC + 2
      DELTA = SI(J1) + SI(J2) + SI(J3) - SI(J4)
      LL=0
      X = ONE
      IF (MAXCHI.LE.0) GOTO 2
      A1 = DBLE(IBCAM + MAXCHI + 1)
      A2 = DBLE(IACBM + MAXCHI + 1)
      A3 = DBLE(MINCHI + MAXCHI + 1)
      B1 = DBLE(IABMC - MAXCHI)
      B2 = DBLE(IAAM - MAXCHI)
      B3 = DBLE(IBBM - MAXCHI)
      DO 1 ICHI = 1, MAXCHI
      XCHI = DBLE(ICHI)
      XB = (A1 - XCHI) * (A2 - XCHI) * (A3 - XCHI)
1     X = ONE - (B1 + XCHI) * (B2 + XCHI) * (B3 + XCHI) * X / XB
      IF(X) 4,3,2
4     X = -X
      LL = 1
2     X = DLOG(X) - SI(IABMC+1) - SI(IAAM+1) - SI(IBBM+1) - SI(IBCAM+1)
     :            - SI(IACBM+1) - SI(MINCHI+1)
      X = TWO*X + SI(IAPAM) + SI(IAMAM) + SI(IBPBM) + SI(IBMBM) +
     :           SI(ICPCM) + SI(ICMCM) + DELTA
      X = DEXP( X / TWO)
      L = LL + MINCHI + NINT(B - A + CM)
      IF( 2 * (L/2) .NE. L) X = -X
3     XF3J = X
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      function plgndr(l,m,x)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause 'bad arguments'
      pmm=1.
      if(m.gt.0) then
         somx2=sqrt((1.-x)*(1.+x))
         fact=1.
         do 11 i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
 11      continue
      endif
      if(l.eq.m) then
         plgndr=pmm
      else
         pmmp1=x*(2*m+1)*pmm
         if(l.eq.m+1) then
            plgndr=pmmp1
         else
            do 12 ll=m+2,l
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
 12         continue
            plgndr=pll
         endif
      endif
      return
      end

C-------------------------------------------------------------------

      SUBROUTINE YTOK(WVEC,L,N,NOPEN,Y,YK,RUP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ROUTINE TO OBTAIN THE K MATRIX FROM 
C     THE LOG DERIVATIVE MATRIX
C     ON ENTRY,  Y HOLDS THE LOG DERIVATIVE MATRIX (N,N)
C     ON EXIT,   YK HOLDS THE AUGMENTED   K MATRIX (NOPEN,NOPEN)
C
C     THESE ARE INPUTS
      DIMENSION WVEC(N), L(N)     
      DIMENSION Y(N,N),YK(NOPEN,NOPEN)
C
C     WVEC(I)=DQSRT(ABS(ETOT-EINT(I))) 
C     L(I)= ORBITAL QUANTUM NUMBER OF CHANNEL I
C
C     THESE ARE INTERNAL VARIABLES
      DIMENSION ABES(N,N), BBES(N,N), 
     *    SJ(N), SJP(N), SN(N), SNP(N)

C     NOTE: Y IS A TWO DIMENSIONAL ARRAY IN DAPROP
      IF(NOPEN.EQ.0) RETURN
      DO I = 1,NOPEN
         DW = WVEC(I)
         DARG = DW*RUP
         CALL BESOPEN(DARG, L(I), RJ, RN, RJP, RNP)
         ROOTDW = SQRT(DW)
         SJ(I) = RJ/ROOTDW
         SJP(I) = RJP*ROOTDW
         SN(I) = RN/ROOTDW
         SNP(I) = RNP*ROOTDW
      ENDDO 
C 
C
      NCLOSE = N - NOPEN
      IF(NCLOSE.GT.0) THEN
       DO I = (NOPEN+1),N
         DW = WVEC(I)
         DARG = DW*RUP
         CALL BESCLOSED(DARG, L(I), RJ, RN, RJP, RNP)
         ROOTDW = SQRT(DW)
         SJ(I) = RJ/ROOTDW
         SJP(I) = RJP*ROOTDW
         SN(I) = RN/ROOTDW
         SNP(I) = RNP*ROOTDW
       ENDDO 
      ENDIF
C
C     FILL BBES=J'-Y*J
C
      DO J=1,N
        DO I=1,N
          BBES(I,J)=-Y(I,J)*SJ(J)
        ENDDO
        BBES(J,J)=BBES(J,J)+SJP(J)
      ENDDO  
C
C     FILL ABES=N'-Y*N
C
      DO J=1,N
        DO I=1,N
          ABES(I,J)=-Y(I,J)*SN(J)
        ENDDO
        ABES(J,J)=ABES(J,J)+SNP(J)
      ENDDO 
C
C     SOLVE ABES*K=BBES TO GET THE AUGMENTED K-MATRIX
C     ON EXIT BBES CONTAIN THE SOLUTION MATRIX
C
      CALL DGESV(N,N,ABES,N,SNP,BBES,N,IERR)
      IF(IERR.NE.0) THEN
       WRITE (6,901) IERR
       STOP
      ENDIF
      DO J=1,NOPEN
       DO I=1,NOPEN
         YK(I,J)=BBES(I,J)
       ENDDO
      ENDDO
      RETURN    
 901  FORMAT('0***** ERROR IN LINEAR EQUATION SOLVER IN YTOK.',
     1     '  IER =',I4,'.  RUN HALTED.')
      END  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE YTOKAUG(WVEC,L,N,NOPEN,Y,YK,RUP,BESZERO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
C     ROUTINE TO OBTAIN THE K MATRIX FROM 
C     THE LOG DERIVATIVE MATRIX
C     WE GET THE AUGMENTED K-MATRIX
C     IN WHICH THE CC,CO AND OC BLOCKS ARE "SCALED"
C     ACCORDINGLY TO THE DEFINITIONS OF THE
C     SCALED MODIFIFIED BESSEL FUNCTIONS
C     IF AUG.EQ.FALSE WE GET THE OPEN-OPEN K MATRIX
C     ON ENTRY,  Y HOLDS THE LOG DERIVATIVE MATRIX (N,N)
C     ON EXIT,   YK HOLDS THE AUGMENTED   K MATRIX (N,N)
C     
      
      LOGICAL BESZERO
      
C     THESE ARE INPUTS
      DIMENSION WVEC(N), L(N)     
      DIMENSION Y(N,N)
      
C     YK IS USED INTERNALLY AND ON OUPUT CONTAINS THE K-MATRIX
      DIMENSION YK(N,N) 
C     
C     WVEC(I)=DQSRT(ABS(ETOT-EINT(I))) 
C     L(I)= ORBITAL QUANTUM NUMBER OF CHANNEL I
C     
C     THESE ARE INTERNAL VARIABLES
      DIMENSION ABES(N,N), SJ(N), SJP(N), SN(N), SNP(N)
      
      FACTOR=DSQRT(2.0D0/DACOS(-1.0D0))

      IF(NOPEN.EQ.0) RETURN
      NCLOSE = N - NOPEN      
C
      IF(BESZERO) THEN
         DO I = 1,NOPEN
            DW = WVEC(I)
            DARG = DW*RUP
            ROOTDW = SQRT(DW)
            RS=DSIN(DARG)
            RC=DCOS(DARG)
            SJ(I) = RS/ROOTDW
            SJP(I) = RC*ROOTDW
            SN(I) = -RC/ROOTDW
            SNP(I) = RS*ROOTDW
         ENDDO 
         IF(NCLOSE.GT.0) THEN
            LZERO=0
            DO I = (NOPEN+1),N
               DW = WVEC(I)
               DARG = DW*RUP
               RJ=FACTOR*(1.0D0-DEXP(-2.0D0*DARG))/2.0D0
               RN=-1.0D0/FACTOR
               RJP=FACTOR*(1.0D0+DEXP(-2.0D0*DARG))/2.0D0
               RNP=1.0D0/FACTOR
C               write(*,*)'IN YTOKAUG'
C               write(*,*) RJ,RN,RJP,RNP
C               CALL BESCLOSED(DARG, LZERO, RJ, RN, RJP, RNP)
C               write(*,*) RJ,RN,RJP,RNP 
               ROOTDW = SQRT(DW)
               SJ(I) = RJ/ROOTDW
               SJP(I) = RJP*ROOTDW
               SN(I) = RN/ROOTDW
               SNP(I) = RNP*ROOTDW
            ENDDO
         ENDIF
C     
      ELSE
C     
         DO I = 1,NOPEN
            DW = WVEC(I)
            DARG = DW*RUP
            CALL BESOPEN(DARG, L(I), RJ, RN, RJP, RNP)
            ROOTDW = SQRT(DW)
            SJ(I) = RJ/ROOTDW
            SJP(I) = RJP*ROOTDW
            SN(I) = RN/ROOTDW
            SNP(I) = RNP*ROOTDW
         ENDDO     
         IF(NCLOSE.GT.0) THEN
            DO I = (NOPEN+1),N
               DW = WVEC(I)
               DARG = DW*RUP
               CALL BESCLOSED(DARG, L(I), RJ, RN, RJP, RNP)
               ROOTDW = SQRT(DW)
               SJ(I) = RJ/ROOTDW
               SJP(I) = RJP*ROOTDW
               SN(I) = RN/ROOTDW
               SNP(I) = RNP*ROOTDW
            ENDDO 
         ENDIF
      ENDIF
C     
C     
C     
      CALL DSYFIL('U',N,Y,N)
C     
C     FILL BBES=J'-Y*J. WE NAMED IT YK
C     
      DO J=1,N
         DO I=1,N
            YK(I,J)=-Y(I,J)*SJ(J)
         ENDDO
         YK(J,J)=YK(J,J)+SJP(J)
      ENDDO  
C     
C     FILL ABES=N'-Y*N
C     
      DO J=1,N
         DO I=1,N
            ABES(I,J)=-Y(I,J)*SN(J)
         ENDDO
         ABES(J,J)=ABES(J,J)+SNP(J)
      ENDDO 
C     
C     SOLVE ABES*K=BBES TO GET THE AUGMENTED K-MATRIX
C     ON EXIT BBES CONTAIN THE SOLUTION MATRIX
C     
      CALL DGESV(N,N,ABES,N,SNP,YK,N,IER)
      IF(IER.NE.0) THEN
         WRITE (6,901) IER
         STOP
      ENDIF
      RETURN    
 901  FORMAT('0***** ERROR IN LINEAR EQUATION SOLVER IN YTOK.',
     1     '  IER =',I4,'.  RUN HALTED.')
      END  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Subroutine besopen(x,n,rj,ry,rjp,ryp)     
C
C-------------------------------------------------------------
C     
C     THIS SUBROUTINES COMPUTE ODINARY AND MODIFIED 
C     RICCATI-BESSEL AND RICCATI-NEUMANN 
C     FUNCTIONS AND THEIR DERIVATIVES.
C 
C     - OPEN CHANNELS:
C       RICCATI-BESSEL (RJ) AND RICCATI-NEUMANN (RY)
C       ARE DEFINED ACCORDINGLY TO ABRAMOWITZ&STEGUN
C            F_N(Z)=Z*FSPH_N(Z) 
C       WHERE 
C            FSPH_N(Z)=DSQRT(PI/2/Z)*F_(N+0.5)(Z)
C       ARE THE USUAL SPHERICAL FUNCTIONS. 
C       THE WRONSKIAN W{RJ,RY} IS 1.
C
C     - CLOSED CHANNELS:
C       SCALED RICCATI-BESSEL (RI) AND RICCATI-NEUMANN (RK)
C       ARE DEFINED AS
C       RI_N(Z) =  DSQRT(Z)*I_(N+0.5) (Z) * EXP(-Z)
C       RK_N(Z) = -DSQRT(Z)*K_(N+0.5) (Z) * EXP(Z)
C       RIP_N(Z)= (RI_N(Z) * EXP(Z))' * EXP(-Z)
C       RKP_N(Z)= (RK_N(Z) * EXP(-Z))' * EXP(Z)
C       THE WRONSKIAN W{RI*EXP(Z),RK*EXP(-Z)} 
C       (NON-SCALED QUANTITIES) IS 1.
C       THE SCALING AVOIDS NUMERICAL UNDERFLOWS AND OVERFLOWS
C
C-----------------------------------------------------------------  
      implicit double precision(a-h,o-z)
      PARAMETER (XMIN=2.,PI=3.141592653589793d0)
c     
c     PI is used for consistency with bessjy
c     
      xnu=dble(n)+0.50d0
      dsx=dsqrt(x)
      factor=dsqrt(pi/2.0d0)
      call bessjy(x,xnu,rj,ry,rjp,ryp)
      rjp=factor*(1.0d0/2.0d0/dsx*rj+dsx*rjp)
      ryp=factor*(1.0d0/2.0d0/dsx*ry+dsx*ryp)
      rj=factor*dsx*rj
      ry=factor*dsx*ry 
      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Subroutine besclosed(x,n,ri,rk,rip,rkp)
      implicit double precision(a-h,o-z)
      PARAMETER (XMIN=2.,PI=3.141592653589793d0)
c
c     change XMIN accordingly to bessik
c     PI is used for consistency with bessik
c
      xnu=dble(n)+0.50d0
      dsx=dsqrt(x)
      call bessik(x,xnu,ri,rk,rip,rkp)
      rip=(1.0d0/2.0d0/dsx*ri+dsx*rip)
      rkp=-(1.0d0/2.0d0/dsx*rk+dsx*rkp)
      ri=dsx*ri
      rk=-dsx*rk
c
c     bessik computes scaled quantities only for x.ge.XMIN
c
      if(x.lt.XMIN) then
       scale=dexp(x)
       ri=ri/scale
       rip=rip/scale
       rk=rk*scale
       rkp=rkp*scale
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp) 
      INTEGER MAXIT 
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI 
      PARAMETER(EPS=1.e-16,FPMIN=1.e-30,MAXIT=1000000,XMIN=2., 
     *     PI=3.141592653589793d0) 

C     USES beschb
C     Return the Bessel functions rj=J_v, ry=Y_v and their derivatives
C     rjp=J'_v, ryp=Y'_v for positive x and for xvu=v.ge.0.
C     The relative accuracy is whithin one or two digits of EPS, except
C     near a zero of one of the functions, where EPS controls its absolute 
C     accuracy. FPMIN is a number close to the machine'smallest floating number.
C     All internal arithmetic is in souble precision. To convert 
C     the entire routine
C     to double precision change REAL declaration above and 
C     decrease EPS to 10**-16.
C     Also convert beschb.

      INTEGER i,isign,l,nl 
      DOUBLE PRECISIONa,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli, 
     *     dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,
     *     p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1, 
     *     rymu,rymup,rytemp,sum,sum1,temp,w,x2,xi,xi2,xmu,xmu2
      
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      
      if(x.lt.XMIN)then 
         nl=int(xnu+.5d0) 
      else
         nl=max(0,int(xnu-x+1.5d0)) 
      endif 
         xmu=xnu-nl 
         xmu2=xmu*xmu 
         xi=1.d0/x
         xi2=2.d0*xi
         w=xi2/PI
         isign=1 
         h=xnu*xi
         if(h.lt.FPMIN)h=FPMIN 
         b=xi2*xnu 
         d=0.d0 
         c=h 
         do i=1,MAXIT
            b=b+xi2 
            d=b-d 
            if(abs(d).lt.FPMIN)d=FPMIN 
            c=b-1.d0/c
            if(abs(c).lt.FPMIN)c=FPMIN 
            d=1.d0/d 
            del=c*d 
            h=del*h
            if(d.lt.0.d0)isign=-isign 
            if(abs(del-1.d0).lt.EPS)goto 1 
         enddo 

c         pause 'x too large in bessjy; try asymptotic expansion' 
         write(6,*)'x too large in bessjy; trying asymptotic expansion' 
         wcos=dcos(x-pi*.5d0*xnu-pi*.25d0)
         wsin=dsin(x-pi*.5d0*xnu-pi*.25d0)
         wk=sqrt(2.d0/(x*pi))
         wkk=x*sqrt(2.d0*pi*x)
         rj=wk*wcos
         ry=wk*wsin
         rjp=(-1.d0/wkk)*wcos-ry
         ryp=(-1.d0/wkk)*wsin+rj
         goto 4
1        continue
         
         rjl=isign*FPMIN 
         rjpl=h*rjl 
         rjl1=rjl 
         rjp1=rjpl 
         fact=xnu*xi 
         do l=nl,1,-1
            rjtemp=fact*rjl+rjpl 
            fact=fact-xi 
            rjpl=fact*rjtemp-rjl 
            rjl=rjtemp
         enddo  
         if(rjl.eq.0.d0)rjl=EPS 
         f=rjpl/rjl
         if(x.lt.XMIN) then
            x2=.5d0*x 
            pimu=PI*xmu 
            if(abs(pimu).lt.EPS)then
               fact=1.d0 
            else
               fact=pimu/sin(pimu) 
            endif 
            d=-log(x2) 
            e=xmu*d 
            if(abs(e).lt.EPS)then
               fact2=1.d0 
            else
               fact2=sinh(e)/e 
            enDif 
            call beschb(xmu,gam1,gam2,gampl,gammi)  
            ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d) 
            e=exp(e)
            p=e/(gampl*PI) 
            q=1.d0/(e*PI*gammi) 
            pimu2=0.5d0*pimu
            if(abs(pimu2).lt.EPS)then
               fact3=1.d0 
            else
               fact3=sin(pimu2)/pimu2 
            endif 
            r=PI*pimu2*fact3*fact3 
            c=1.d0 
            d=-x2*x2
            sum=ff+r*q 
            sum1=p
            do  i=1,MAXIT
               ff=(i*ff+p+q)/(i*i-xmu2) 
               c=c*d/i 
               p=p/(i-xmu) 
               q=q/(i+xmu)
               del=c*(ff+r*q) 
               sum=sum+del 
               del1=c*p-i*del 
               sum1=sum1+del1
               if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2 
            enddo  
            pause 'bessy series failed to converge'
2           continue
            rymu=-sum 
            ry1=-sum1*xi2 
            rymup=xmu*xi*rymu-ry1 
            rjmu=w/(rymup-f*rymu)
         else 
            a=.25d0-xmu2 
            p=-.5d0*xi 
            q=1.d0 
            br=2.d0*x 
            bi=2.d0
            fact=a*xi/(p*p+q*q) 
            cr=br+q*fact 
            ci=bi+p*fact 
            den=br*br+bi*bi
            dr=br/den 
            di=-bi/den 
            dlr=cr*dr-ci*di 
            dli=cr*di+ci*dr 
            temp=p*dlr-q*dli
            q=p*dli+q*dlr 
            p=temp 
            do i=2,MAXIT
               a=a+2*(i-1) 
               bi=bi+2.d0 
               dr=a*dr+br 
               di=a*di+bi
               if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN 
               fact=a/(cr*cr+ci*ci)
               cr=br+cr*fact 
               ci=bi-ci*fact 
               if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
               den=dr*dr+di*di 
               dr=dr/den 
               di=-di/den 
               dlr=cr*dr-ci*di 
               dli=cr*di+ci*dr
               temp=p*dlr-q*dli 
               q=p*dli+q*dlr 
               p=temp
               if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3 
            enddo  
            pause 'cf2 failed in bessjy' 
3            continue
            gam=(p-f)/q 
            rjmu=sqrt(w/((p-f)*gam+q))
            rjmu=sign(rjmu,rjl) 
            rymu=rjmu*gam 
            rymup=rymu*(p+q/gam)
            ry1=xmu*xi*rymu-rymup 
         endif 
         fact=rjmu/rjl
         rj=rjl1*fact  
         rjp=rjp1*fact 
         do  i=1,nl
            rytemp=(xmu+i)*xi2*ry1-rymu 
            rymu=ry1 
            ry1=rytemp 
         enddo  
         ry=rymu
         ryp=xnu*xi*rymu-ry1 
4        return 
         END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp) 
      INTEGER MAXIT 
      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN 
      DOUBLE PRECISION EPS,FPMIN,PI 
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=1000000,XMIN=2., 
     *     PI=3.141592653589793d0) 
      
C     USES beschb
C     Returns the modified Bessel functions ri=I_v, rk=K_v and their derivatives
C     rip=I'_v, rkp=K'_v for positive x and for xnu=v.ge.0. The relative accuracy
C     is whithin one or two significant digits of EPS. FPMIN is a number close to 
C     the machine's smallest floating point number. All internal arithmetic is in 
C     double precision. To convert the entire routine to double precision, change
C     REAL declaration above and decrease EPS to 10**-16.
C     Also convert the subroutine beschb. 
      
      INTEGER i,l,nl 
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact, 
     *     fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2, 
     *     qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup, 
     *     rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl 
      xmu2=xmu*xmu 
      xi=1.d0/x 
      xi2=2.d0*xi 
      h=xnu*xi 
      if(h.lt.FPMIN)h=FPMIN 
      b=xi2*xnu 
      d=0.d0 
      c=h 
      do i=1,MAXIT
         b=b+xi2 
         d=1.d0/(b+d) 
         c=b+1.d0/c 
         del=c*d 
         h=del*h
         if(abs(del-1.d0).lt.EPS)goto 1 
      enddo  

      pause 'x too large in bessik;try asymptotic expansion' 
 1    continue
      ril=FPMIN 
      ripl=h*ril
      ril1=ril
      rip1=ripl 
      fact=xnu*xi 
      do l=nl,1,-1
         ritemp=fact*ril+ripl 
         fact=fact-xi 
         ripl=fact*ritemp+ril 
         ril=ritemp
      enddo  
      f=ripl/ril  
      if(x.lt.XMIN)then 
         x2=.5d0*x 
         pimu=PI*xmu 
         if(abs(pimu).lt.EPS)then
            fact=1.d0 
         else
            fact=pimu/sin(pimu) 
         endif 
         d=-log(x2) 
         e=xmu*d 
         if(abs(e).lt.EPS)then
            fact2=1.d0 
         else
            fact2=sinh(e)/e 
         endif 
         call beschb(xmu,gam1,gam2,gampl,gammi)  
         ff=fact*(gam1*cosh(e)+gam2*fact2*d) 
         sum=ff 
         e=exp(e)
         p=0.5d0*e/gampl 
         q=0.5d0/(e*gammi)  
         c=1.d0 
         d=x2*x2 
         sum1=p 
         do i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2) 
            c=c*d/i 
            p=p/(i-xmu) 
            q=q/(i+xmu) 
            del=c*ff
            sum=sum+del 
            del1=c*(p-i*ff)
            sum1=sum1+del1 
            if(abs(del).lt.abs(sum)*EPS)goto 2 
         enddo  
         pause 'bessk series failed to converge' 
 2       continue
         
         rkmu=sum 
         rk1=sum1*xi2 
      else  
         b=2.d0*(1.d0+x)
         d=1.d0/b 
         delh=d 
         h=delh 
         q1=0.d0  
         q2=1.d0 
         a1=.25d0-xmu2 
         c=a1 
         q=c 
         a=-a1 
         s=1.d0+q*delh 
         do i=2,MAXIT
            a=a-2*(i-1) 
            c=-a*c/i 
            qnew=(q1-b*q2)/a 
            q1=q2 
            q2=qnew 
            q=q+c*qnew
            b=b+2.d0 
            d=1.d0/(b+a*d) 
            delh=(b*d-1.d0)*delh 
            h=h+delh 
            dels=q*delh
            s=s+dels 
            if(abs(dels/s).lt.EPS)goto 3 
         enddo  
         pause 'bessik: failure to converge in cf2' 
 3       continue
         h=a1*h 
c
c     To compute scaled quantities
c
c         rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s 
c
         rkmu=sqrt(PI/(2.d0*x))/s 
         rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif 
      rkmup=xmu*xi*rkmu-rk1 
      rimu=xi/(f*rkmu-rkmup)  
      ri=(rimu*ril1)/ril  
      rip=(rimu*rip1)/ril 
      do i=1,nl 
         rktemp=(xmu+i)*xi2*rk1+rkmu 
         rkmu=rk1 
         rk1=rktemp 
      enddo  
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1 
      return 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi) 
      INTEGER NUSE1,NUSE2 
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x 
      PARAMETER (NUSE1=7,NUSE2=8) 

C     USES chebev
C     Evluata Gamma1 and Gamma2 by Chebischev expansion for abs(x).le.0.5.
C     Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
C     If converting to double precision, set NUSE1=7, NUSE2=8.

      DOUBLE PRECISION xx,c1(7),c2(8),chebev 
      SAVE c1,c2 
      DATA c1/-1.142022680371168d0,6.5165112670737d-3, 
     *     3.087090173086d-4,-3.4706269649d-6,6.9437664d-9, 
     *     3.67795d-11,-1.356d-13/
      
      DATA c2/1.843740587300905d0,-7.68528408447867d-2, 
     *     1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8, 
     *     2.423096d-10,-1.702d-13,-1.49d-15/
      
      xx=8.d0*x*x-1.d0 
      gam1=chebev(-1.0d0,1.0d0,c1,NUSE1,xx)
      gam2=chebev(-1.0d0,1.0d0,c2,NUSE2,xx) 
      gampl=gam2-x*gam1 
      gammi=gam2+x*gam1
      return 
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function chebev(a,b,c,m,x)
      implicit double precision (a-h,o-z)
      dimension c(m)
      if ((x-a)*(x-b).gt.0.d0) pause 'x not in range.'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
 11   continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE KTOY(WVEC,L,N,NOPEN,YK,Y,RUP,BESZERO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     THESE ARE INPUTS
      LOGICAL BESZERO
      DIMENSION WVEC(N), L(N)     
      DIMENSION Y(N,N)

      DIMENSION YK(N,N) 
C
C     WVEC(I)=DQSRT(ABS(ETOT-EINT(I))) 
C     L(I)= ORBITAL QUANTUM NUMBER OF CHANNEL I
C
C     THESE ARE INTERNAL VARIABLES
      DIMENSION ABES(N,N), SJ(N), SJP(N), SN(N), SNP(N)

      FACTOR=DSQRT(2.0D0/DACOS(-1.0D0))
      IF(NOPEN.EQ.0) RETURN
      NCLOSE = N - NOPEN
C
      IF(BESZERO) THEN
       DO I = 1,NOPEN
         DW = WVEC(I)
         DARG = DW*RUP
         ROOTDW = SQRT(DW)
         RS=DSIN(DARG)
         RC=DCOS(DARG)
         SJ(I) = RS/ROOTDW
         SJP(I) = RC*ROOTDW
         SN(I) = -RC/ROOTDW
         SNP(I) = RS*ROOTDW
       ENDDO
       IF(NCLOSE.GT.0) THEN
          LZERO=0
          DO I = (NOPEN+1),N
             DW = WVEC(I)
             DARG = DW*RUP
             RJ=FACTOR*(1.0D0-DEXP(-2.0D0*DARG))/2.0D0
             RN=-1.0D0/FACTOR
             RJP=FACTOR*(1.0D0+DEXP(-2.0D0*DARG))/2.0D0
             RNP=1.0D0/FACTOR
C             CALL BESCLOSED(DARG, LZERO, RJ, RN, RJP, RNP)
             ROOTDW = SQRT(DW)
             SJ(I) = RJ/ROOTDW
             SJP(I) = RJP*ROOTDW
             SN(I) = RN/ROOTDW
             SNP(I) = RNP*ROOTDW
          ENDDO 
       ENDIF
      ELSE
         DO I = 1,NOPEN
            DW = WVEC(I)
            DARG = DW*RUP
            CALL BESOPEN(DARG, L(I), RJ, RN, RJP, RNP)
            ROOTDW = SQRT(DW)
            SJ(I) = RJ/ROOTDW
            SJP(I) = RJP*ROOTDW
            SN(I) = RN/ROOTDW
            SNP(I) = RNP*ROOTDW
         ENDDO 
         IF(NCLOSE.GT.0) THEN
            DO I = (NOPEN+1),N
               DW = WVEC(I)
               DARG = DW*RUP
               CALL BESCLOSED(DARG, L(I), RJ, RN, RJP, RNP)
               ROOTDW = SQRT(DW)
               SJ(I) = RJ/ROOTDW
               SJP(I) = RJP*ROOTDW
               SN(I) = RN/ROOTDW
               SNP(I) = RNP*ROOTDW
            ENDDO 
         ENDIF
      ENDIF
C
C     
C     CALL DSYFIL('U',N,YK,N)
C     
C     FILL ABES=PSI_TRANS=J-K*N
C     
      DO J=1,N
         DO I=1,N
            ABES(I,J)=-YK(I,J)*SN(J)
         ENDDO
         ABES(J,J)=ABES(J,J)+SJ(J)
      ENDDO  
C     
C     FILL PSI_TRANS'=J'-K*N'
C     
      DO J=1,N
         DO I=1,N
            Y(I,J)=-YK(I,J)*SNP(J)
         ENDDO
         Y(J,J)=Y(J,J)+SJP(J)
      ENDDO 
C     
C     SOLVE (J-K*N)*YTRAN=(J'-K*N') TO GET THE Y-MATRIX
C     ON EXIT Y  CONTAIN THE SOLUTION MATRIX
C     
      CALL DGESV(N,N,ABES,N,SNP,Y,N,IER)
      IF(IER.NE.0) THEN
         WRITE (6,901) IER
         STOP
      ENDIF
      RETURN    
 901  FORMAT('0***** ERROR IN LINEAR EQUATION SOLVER IN YTOK.',
     1     '  IER =',I4,'.  RUN HALTED.')
      END  
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      SUBROUTINE KTOK(WVEC,L,NINIT,N,NOPEN,YK,RUP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      LOGICAL BESZERO
      DIMENSION WVEC(N), L(N)
      DIMENSION YK(NINIT,*)
      DIMENSION YK1(N,N),Y(N,N)
      
      DO J=1,N
         DO I=1,N
            YK1(I,J)=YK(I,J)
         ENDDO
      ENDDO
      
C     THE BESZERO-K MATRIX IS CONVERTED INTO Y MATRIX
      BESZERO=.TRUE.
      CALL KTOY(WVEC,L,N,NOPEN,YK1,Y,RUP,BESZERO)
      
C     THE Y MATRIX IS CONVERTED INTO USUAL K MATRIX
      BESZERO=.FALSE.
      CALL YTOKAUG(WVEC,L,N,NOPEN,Y,YK1,RUP,BESZERO)
      
      DO J=1,N
         DO I=1,N
            YK(I,J)=YK1(I,J)
         ENDDO
      ENDDO
      
C      WRITE(*,*) '                         '
C      WRITE(*,*) 'CONVERSION TO USUAL K MATRIX AT R=',RUP,' A.U.'
C      WRITE(*,*) '                         '
      RETURN 
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ktos(xk,pr,pi,n)
      implicit double precision(a-h,o-z) 
c Calculation of the scattering matrix S: S:=(1+iK)*inv(1-iK).
c K:= reactance matrix, real and symmetrical.
c S is complex.
c Check that:
c real part of S: R:=(1+K**2)*inv(1+K**2),
c imaginary part of S: I:=2*K*inv(1+K**2). 
c xk is K, xk2 is K**2,
c umask2 is 1+K**2, umenosk2 is 1-K**2,
c pr is R and pi is I.
      dimension xk(n,n), xk2(n,n)
      dimension umask2(n,n), umenosk2(n,n)
      dimension pr(n,n), pi(n,n)

c K**2 matrix

      call dsymm('l','l',n,n,1.0d0,xk,n,xk,n,0.0d0,xk2,n)

      do i=1,n
      do j=i+1,n
      xk2(i,j)=xk2(j,i)
      enddo
      enddo 

c 1+K**2 matrix

      do i=1,n
      do j=1,n
      if (i.eq.j) then
      umask2(i,j)=1.0d0+xk2(i,j) 
      else
      umask2(i,j)=xk2(i,j)
      endif
      enddo
      enddo

c inverse of 1+K**2

      call syminv(umask2,n,n)

      do i=1,n
      do j=i+1,n
      umask2(i,j)=umask2(j,i)
      enddo
      enddo

c 1-K*K matrix

      do i=1,n
      do j=1,n
      if(i.eq.j)then
      umenosk2(i,j)=1.0d0-xk2(i,j)
      else
      umenosk2(i,j)=-xk2(i,j)
      endif
      enddo
      enddo

c real part of S, R=(1-K**2)*inv(1+K**2) 

      call dsymm('l','l',n,n,1.0d0,umenosk2,n,umask2,n,0.0d0,pr,n)

      do i=1,n
      do j=i+1,n
      pr(i,j)=pr(j,i)
      enddo
      enddo

c imaginary part of S, I=2*K*inv(1+K**2)

      call dsymm('l','l',n,n,2.0d0,xk,n,umask2,n,0.0d0,pi,n)

      do i=1,n
      do j=i+1,n
      pi(i,j)=pi(j,i)
      enddo
      enddo
      return
      end








