!!!*************************************************************
! 文件/File: spl1d1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: spl1d1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  SPL1D1
C***********************************************************************
C
      SUBROUTINE SPL1D1 (N,X,F,W,IOP,IJ,A,B,C)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C     WHERE N = NUMBER OF POINTS IN THE INTERPOLATION
C           X = ORIGIN OF TABLE OF INDEPENDENT VARIALE
C           F = ORIGIN OF TABLE OF DEPENDENT VARIABLE
C           W = AN ARRAY OF DIMENSION N WHICH CONTAINS THE CALCULATED
C               SECOND DERIVATIVES UPON RETURN
C         IOP = AN ARRAY OF DIMENSION 2 WHICH CONTAINS COMBINATIONS OF
C               THE INTEGERS 1 THRU 5 USED TO SPECIFY THE BOUNDARY
C               CONDITIONS
C     IF IOP(1)=1 AND IOP(2)=1, W(1) AND W(K) ARE THE VALUES OF THE
C     SECOND DERIVATIVE AT X(1) AND X(K), RESPECTIVELY.
C     IF IOP(1)=2 AND IOP(2)=2, W(1) DETERMINES F""(X(1)) BY THE
C     RELATION F""(X(1))=W(1)*F""(X(2)), AND W(K) DOES F""(X(K)) BY THE
C     RELATION F""(X(N))=W(K)*F""(X(N-1)).
C     IF IOP(1)=3 AND IOP(2)=3, W(1) AND W(K) ARE THE VALUES OF THE
C     FIRST DERIVATIVE AT X(1) AND X(K), RESPECTIVELY.
C     IF IOP(1)=4 AND IOP(2)=4, THE RERIODIC BOUNDARY CONDITIONS ARE
C     ALLOWED, F""(1)=F""(N), F(1)=F(N).
C     IF IOP(1)=5 AND IOP(2)=5, THE FIRST DERIVATIVES AT X(1) AND X(N)
C     ARE CALCULATED BY USING A DIFFERENTIATED FOUR-POINT LAGRANGIAN
C     INTERPOLATION FORMULA EVALUATED AT X(1) AND X(N), RESPECTIVELYR.
C          IJ = SPACING IN THE F AND W TABLES
C       A,B,C = ARRAYS OF DIMENSION N USED FOR TEMPORARY STORAGE
C
      DIMENSION IOP(2),X(N),F(N),W(N),A(N),B(N),C(N)
      SAVE BOB, BILL, EPS                                               TCA1097
      DATA BOB / 0.0D0 /, BILL / 0.0D0 /, EPS/ 1.0D-11 /
      K = N-1
      SXIV = (1.0D0/6.0D0)-EPS
      THIV = 2.0D0*SXIV
      A(2) = -(X(2)-X(1))*SXIV
      B(2) = (X(3)-X(1))*THIV
      W(IJ+1) = (F(2*IJ+1)-F(IJ+1))/(X(3)-X(2))-(F(IJ+1)-F(1))/(X(2)-X(1
     *   ))
      IF (N-3) 10, 30, 10
   10 DO 20 I = 3, K
         M = (I-1)*IJ+1
         J1 = M+IJ
         J2 = M-IJ
         CON = (X(I+1)-X(I-1))*THIV
         DIFX = X(I)-X(I-1)
         DON = DIFX*SXIV
         BIMI = 1.0D0/B(I-1)
         B(I) = CON-(DON**2)*BIMI
         E = (F(J1)-F(M))/(X(I+1)-X(I))-(F(M)-F(J2))/DIFX
         W(M) = E-(DON*W(J2))*BIMI
         A(I) = -(DON*A(I-1))*BIMI
   20 CONTINUE
   30 K1 = (N-2)*IJ+1
      BNMI = 1.0D0/B(N-1)
      C(N-1) = -(X(N)-X(N-1))*SXIV*BNMI
      W(K1) = W(K1)*BNMI
      A(N-1) = A(N-1)*BNMI
      K2 = K-1
      IF (N-3) 40, 60, 40
   40 DO 50 I = 2, K2
         J = N-I
         CON = (X(J+1)-X(J))*SXIV
         BJI = 1.0D0/B(J)
         A(J) = (A(J)-CON*A(J+1))*BJI
         C(J) = -(CON*C(J+1))*BJI
         K3 = (J-1)*IJ+1
         M = K3+IJ
         W(K3) = (W(K3)-CON*W(M))*BJI
   50 CONTINUE
   60 K4 = (N-1)*IJ+1
      IF (IOP(1)-5) 70, 90, 70
   70 C1 = W(1)
      IF (IOP(2)-5) 80, 110, 80
   80 C2 = W(K4)
      GO TO 130
   90 IF (N-4) 570, 100, 100
  100 A1 = X(1)-X(2)
      A2 = X(1)-X(3)
      A3 = X(1)-X(4)
      A4 = X(2)-X(3)
      A5 = X(2)-X(4)
      A6 = X(3)-X(4)
      W(1) = F(1)*(1.0D0/A1+1.0D0/A2+1.0D0/A3)-A2*A3*F(IJ+1)/(A1*A4*A5)+
     *   A1*A3*F(2*IJ+1)/(A2*A4*A6)-A1*A2*F(3*IJ+1)/(A3*A5*A6)
      GO TO 70
  110 IF (N-4) 570, 120, 120
  120 B1 = X(N)-X(N-3)
      B2 = X(N)-X(N-2)
      B3 = X(N)-X(N-1)
      B4 = X(N-1)-X(N-3)
      B5 = X(N-1)-X(N-2)
      B6 = X(N-2)-X(N-3)
      L1 = K4-IJ
      L2 = L1-IJ
      L3 = L2-IJ
      W(K4) = -B2*B3*F(L3)/(B6*B4*B1)+B1*B3*F(L2)/(B6*B5*B2)-B1*B2*F(L1)
     *   /(B4*B5*B3)+F(K4)*(1.0D0/B1+1.0D0/B2+1.0D0/B3)
      GO TO 80
  130 DO 160 I = 1, K
         M = (I-1)*IJ+1
  170 MK = IOP(1)
      GO TO (180,210,260,310,260), MK
  180 IF (I-1) 200, 190, 200
  190 A(1) = -1.0D0
      C(1) = 0.0D0
      GO TO 340
  200 BOB = 0.0D0
      GO TO 340
  210 IF (I-1) 230, 220, 230
  220 A(1) = -1.0D0
      C(1) = 0.0D0
      W(1) = 0.0D0
      GO TO 340
  230 IF (I-2) 240, 240, 250
  240 BOB = -C1
      GO TO 340
  250 BOB = 0.0D0
      GO TO 340
  260 IF (I-1) 280, 270, 280
  270 XDTO = X(2)-X(1)
      A(1) = -XDTO*THIV
      C(1) = 0.0D0
      W(1) = -C1+(F(IJ+1)-F(1))/XDTO
      GO TO 340
  280 IF (I-2) 290, 290, 300
  290 BOB = (X(2)-X(1))*SXIV
      GO TO 340
  300 BOB = 0.0D0
      GO TO 340
  310 IF (I-1) 330, 320, 330
  320 A(1) = -1.0D0
      C(1) = 1.0D0
      W(1) = 0.0D0
      GO TO 340
  330 BOB = 0.0D0
  340 ML = IOP(2)
      GO TO (350,380,430,480,430), ML
  350 IF (I-1) 370, 360, 370
  360 A(N) = 0.0D0
      C(N) = -1.0D0
      GO TO 140
  370 BILL = 0.0D0
      GO TO 140
  380 IF (I-1) 400, 390, 400
  390 A(N) = 0.0D0
      C(N) = -1.0D0
      W(K4) = 0.0D0
      GO TO 140
  400 IF (I-K) 420, 410, 420
  410 BILL = -C2
      GO TO 140
  420 BILL = 0.0D0
      GO TO 140
  430 IF (I-1) 450, 440, 450
  440 A(N) = 0.0D0
      C(N) = (X(N-1)-X(N))*THIV
      W(K4) = C2-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 140
  450 IF (I-K) 470, 460, 470
  460 BILL = (X(N)-X(N-1))*SXIV
      GO TO 140
  470 BILL = 0.0D0
      GO TO 140
  480 IF (I-1) 500, 490, 500
  490 A(N) = 0.0D0
      C(N) = (X(N-1)+X(1)-X(N)-X(2))*THIV
      W(K4) = (F(IJ+1)-F(1))/(X(2)-X(1))-(F(K4)-F(K1))/(X(N)-X(N-1))
      GO TO 140
  500 IF (I-2) 520, 510, 520
  510 BILL = (X(2)-X(1))*SXIV
      GO TO 140
  520 IF (I-K) 540, 530, 540
  530 BILL = (X(N)-X(N-1))*SXIV
      GO TO 140
  540 BILL = 0.0D0
      GO TO 140
  140    IF (I-1) 150, 160, 150
  150    W(1) = W(1)-BOB*W(M)
         W(K4) = W(K4)-BILL*W(M)
         A(1) = A(1)-BOB*A(I)
         A(N) = A(N)-BILL*A(I)
         C(1) = C(1)-BOB*C(I)
         C(N) = C(N)-BILL*C(I)
  160 CONTINUE
  550 CON = A(1)*C(N)-C(1)*A(N)
      D1 = -W(1)
      D2 = -W(K4)
      W(1) = (D1*C(N)-C(1)*D2)/CON
      W(K4) = (A(1)*D2-D1*A(N))/CON
      DO 560 I = 2, K
         M = (I-1)*IJ+1
         W(M) = W(M)+A(I)*W(1)+C(I)*W(K4)
  560 CONTINUE
      GO TO 580
  570 WRITE (6,1000)
  580 RETURN
C
1000  FORMAT(2X,T5,'Warning: In SPL1D1 N is less than 4, the ',
     *             'results are incorrect')
C
      END
