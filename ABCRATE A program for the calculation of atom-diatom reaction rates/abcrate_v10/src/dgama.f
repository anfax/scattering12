!!!*************************************************************
! 文件/File: dgama.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: dgama.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  DGAMA
C***********************************************************************
C
      FUNCTION  DGAMA(ZZ)                                               GCL0992
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C     MODIFIED 4/4/88 BY DWS TO BE ACCURATE FOR ZZ GT 3.                DWS
C  THIS IS A PROCEDURE THAT EVALUATES GAMMA(Z) FOR                      GD.426
C     0 LT Z LE 3 TO 16 SIGNIFICANT FIGURES                             GD.427
C    IT IS BASED ON A CHEBYSHEV-TYPE POLYNOMIAL                         GD.428
C   APPROXIMATION GIVEN IN H. WERNER AND R. COLLINGE, MATH. COMPUT.     GD.429
C    15 (1961), PP. 195-97.                                             GD.430
C   APPROXIMATIONS TO THE GAMMA FUNCTION, ACCURATE UP TO 18 SIGNIFICANT GD.431
C   DIGITS, MAY BE FOUND IN THE PAPER QUOTED ABOVE                      GD.432
C                                                                       GD.433
C                                                                       GD.434
C                                                                       GD.435
C   This function has been renamed DGAMA (previously called DGAMMA)     GCL0992
C   to avoid conflicts with the FORTRAN intrinsic function DGAMMA.      GCL0992
C   All E format has been converted to D format, this change was not    GCL1092
C   initialed in columns 73-80.                                         GCL1092
C
      DIMENSION  A(18)                                                  GD.436
      PREFAC=1.0D0                                                      DWS
      Z=ZZ                                                              DWS
   44 CONTINUE                                                          DWS
      IF(Z.LT.3.0D0)GO TO 45                                            DWS
      Z=Z-1.0D0                                                         DWS
      PREFAC=PREFAC*Z                                                   DWS
      GO TO 44                                                          DWS
   45 CONTINUE                                                          DWS
C                                                                       GD.437
       A(1)=1.0D0                                                       GD.438
       A(2)=.4227843350984678D0                                         GD.439
       A(3)=.4118403304263672D0                                         GD.440
      A(4)=.0815769192502609D0                                          GD.441
      A(5)=.0742490106800904D0                                          GD.442
      A(6)=-.0002669810333484D0                                         GD.443
      A(7)=.0111540360240344D0                                          GD.444
      A(8)=-.0028525821446197D0                                         GD.445
      A(9)=.0021036287024598D0                                          GD.446
      A(10)=-.0009184843690991D0                                        GD.447
      A(11)=.0004874227944768D0                                         GD.448
      A(12)=-.0002347204018919D0                                        GD.449
      A(13)=.0001115339519666D0                                         GD.450
      A(14)=-.0000478747983834D0                                        GD.451
      A(15)=.0000175102727179D0                                         GD.452
      A(16)=-.0000049203750904D0                                        GD.453
      A(17)=.0000009199156407D0                                         GD.454
      A(18)=-.0000000839940496D0                                        GD.455
C                                                                       GD.456
C                                                                       GD.457
C                                                                       GD.458
      IF(Z.LE.1.0D0  ) GO TO 10                                         GD.459
      IF(Z.LE.2.0D0  ) GO TO 20                                         GD.460
      T=Z-2.0D0                                                         GD.461
      GO TO 30                                                          GD.462
10    T=Z                                                               GD.463
      GO TO 30                                                          GD.464
20    T=Z-1.0D0                                                         GD.465
30    P=A(18)                                                           GD.466
      DO 40 K1=1,17                                                     GD.467
      K=18-K1                                                           GD.468
      P=T*P+A(K)                                                        GD.469
40    CONTINUE                                                          GD.470
C                                                                       GD.471
      P=P*PREFAC                                                        DWS
      IF(Z.GT.2.0D0  ) GO TO 50                                         GD.472
      IF(Z.GT.1.0D0  ) GO TO 60                                         GD.473
      DGAMA=P/(Z*(Z+1.0D0  ))                                           GCL0992
      RETURN                                                            GD.475
60    DGAMA=P/Z                                                         GCL0992
      RETURN                                                            GD.477
50    DGAMA=P                                                           GCL0992
      RETURN                                                            GD.479
      END                                                               GD.480
