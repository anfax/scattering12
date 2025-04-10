!!!*************************************************************
! 文件/File: bendfunc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: bendfunc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      THIS FILE INCLUDES SUBROUTINES DEVOTED TO CALCULATION    C
C      OF THE MATRIX ELEMENTS ON THE BENDING BASE FUNCTIONS     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATES MATRIX ELEMENTS OF THE DIPOLE MOMENT OF THE 
C             : TRIATOMIC MOLECULE USING A NUMERICAL INTAGRATION 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE BendMatrixElement 
C
      INTEGER K1,K2,V2_1,V2_2,ISurf,I1,Ipts,IXYZ,i0

      include 'dipolsys.h'
      include 'mu.h'

C
      REAL*8 Phi_V_Phi(0:NumPtsRho),CoefMu(3,NumPtsRho)

      REAL*8 BendFuncAB(NumPtsRho),
     2       BendFuncC(NumPtsRho),SQ2

      REAL*8 DipMomentBend(2,2,UpperV2,UpperV2,UpperK,2,15)

      COMMON /DIPBEND/ DipMomentBend

C
      REAL*8  CosGamma(KMAXAB+1,NumPtsRho),
     1        SinGamma(KMAXAB+1,NumPtsRho),SimpsonIntegral

      INTEGER DeltaK,MinDeltaK(3),NDeltaK,theta

      INTEGER NFileA,NFileB,NFileC,NUNIT,KREC

      DATA MinDeltaK /-1,-1,0/

      CHARACTER*80 FILNAM
C
C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C


      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   BENDING MATRIX ELEMENTS '/)



C**********************************************************************C
C                                                                      C
C       OPEN DIRECT ACCESS FILE FOR STORING BENDING WAVE FUNCTIONS     C
C                                                                      C
C**********************************************************************C
C                                                                      C
      NFileA =  81
      NFileB =  82
      NFileC =  83

      NUNIT=NFileA
      OPEN(UNIT=NUNIT,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NumPtsRho*NRECUN,ERR=1009)
      NUNIT=NFileB
      OPEN(UNIT=NUNIT,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NumPtsRho*NRECUN,ERR=1009)
      NUNIT=NFileC
      OPEN(UNIT=NUNIT,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NumPtsRho*NRECUN,ERR=1009)
C                                                                      C
C**********************************************************************C  


C**********************************************************************C
C                                                                      C
C     READ THE RHO DEPENDENT WAVE FUNCTIONS                            C
C                                                                      C
C**********************************************************************C
C                                                                      C
      CALL ReadBendFunction(NFileA,NFileB,NFileC)
C
C     READ THE RHO DEPENDENT COS AND SIN FUNCTIONS 
C
      CALL ReadCosSinGamma(CosGamma,SinGamma)
C                                                                      C
C**********************************************************************C  

      WRITE(6,3002)
3002  FORMAT(1H0,'                 CALCULATIONS '/)
      SQ2 = DSQRT(2.0D+00)

C**********************************************************************C
C                                                                      C
C     START BENDING MATRIX ELEMENTS CALCULATIONS                       C
C                                                                      C
C**********************************************************************C
C                                                                      C


      DO I0=1,15     
      DO ISurf=1,2       ! --- A,B - electronic states  Left  < |
      DO theta = 0,1                 ! --- theta
      DO K1 = 0,KMAXAB
      DO DeltaK = -theta,theta,2    ! --- DeltaK 
         NDeltaK = INT((DeltaK+2)/3)+1
      DO V2_1=0,NV2Max(ISurf)
      DO V2_2=0,NV2MaxC

       DipMomentBend(theta+1,ISurf,V2_1+1,V2_2+1,K1+1,NDeltaK,I0) = 0.0         

      ENDDO    !   -   end of loop ever V2_2 quantum number 
      ENDDO    !   -   end of loop ever V2_1 quantum number 
      ENDDO    !   -   end of loop ever DeltaK = K2-K1  
      ENDDO    !   -   end of loop ever K1  quantum number 
      ENDDO    !   -   theta 
      ENDDO    !   -   A,B - electronic states  Left  < |
      ENDDO    !   -   end of loop OVER ALL DIFERENT Y_1^i * Y_3^j COMBINATIONS





C-------------------------------------------------------------C
C-------------------------------------------------------------C
C     START LOOP OVER ALL DIFERENT Y_1^i * Y_3^j COMBINATIONS C
C                                                             C
      DO I0=1,15     

C     --- X,Y OR Z DIPOL MOMENT COMPONENT 
      DO  IXYZ =1,3  
          READ(NFileMu) ( CoefMu(IXYZ,Ipts), Ipts = 1,NumPtsRho )
      ENDDO    !   -    X,Y OR Z



      DO  ISurf=1,2       ! --- A,B - electronic states  Left  < |

      IF ( ISurf.EQ.1 ) THEN 
         NUNIT = NFileA
      ELSE 
         NUNIT = NFileB
      ENDIF 

C
C     START LOOP OVER K AND V2 BENDING QUANTUM NUMBERS 
C

      DO theta = 0,1                 ! --- theta
      DO K1 = 0,KMAXAB
      DO DeltaK = -theta,theta,2    ! --- DeltaK 
           K2 = K1 + DeltaK
      IF ( K2.GE.0 .AND. K2.LE.KMAXC ) THEN                ! --- K2 <=KMAXC 
C     ---- NDeltaK  - index for distinguishing DeltaK = -1,0,1 => NdeltaK = 1,1,2
           NDeltaK = INT((DeltaK+2)/3)+1

      DO V2_1=0,NV2Max(ISurf)

      KREC=K1*( NV2Max(ISurf)+1 )+V2_1+1
      READ (NUNIT,REC=KREC)  (BendFuncAB(Ipts),Ipts=1,NumPtsRho)

      DO V2_2=0,NV2MaxC

      KREC=K2*( NV2MaxC+1 )+V2_2+1
      READ (NFileC,REC=KREC)  (BendFuncC(Ipts),Ipts=1,NumPtsRho)



C*****************************************************************C
      IF ( ISurf.EQ.1 .AND. theta.EQ.0 ) THEN 

       Phi_V_Phi(0) =  0  
       DO Ipts = 1,NumPtsRho
        Phi_V_Phi(Ipts) = BendFuncAB(Ipts)*
     1                    CosGamma(K1+1,Ipts)*
     2                    CoefMu(3,Ipts)*
     3                    BendFuncC(Ipts)
       ENDDO 

      ELSEIF ( ISurf.EQ.2 .AND. theta.EQ.0 ) THEN 


       Phi_V_Phi(0) =  0
       DO Ipts = 1,NumPtsRho
        Phi_V_Phi(Ipts) =-BendFuncAB(Ipts)*
     1                    SinGamma(K1+1,Ipts)*
     2                    CoefMu(3,Ipts)*
     3                    BendFuncC(Ipts)
       ENDDO 

      ELSEIF ( ISurf.EQ.1 .AND. theta.EQ.1 )  THEN 

       Phi_V_Phi(0) =  0 
       DO Ipts = 1,NumPtsRho
        Phi_V_Phi(Ipts) = BendFuncAB(Ipts)*SQ2*
     1             ( DeltaK*CosGamma(K1+1,Ipts)*CoefMu(2,Ipts) 
     2                    + SinGamma(K1+1,Ipts)*CoefMu(1,Ipts) )*
     3                    BendFuncC(Ipts)
       ENDDO 

      ELSE  ! --  ( ISurf.EQ.2 .AND. theta.EQ.1 )

       Phi_V_Phi(0) =  0 
       DO Ipts = 1,NumPtsRho
        Phi_V_Phi(Ipts) = BendFuncAB(Ipts)*SQ2*
     1             (-DeltaK*SinGamma(K1+1,Ipts)*CoefMu(2,Ipts) 
     2                    + CosGamma(K1+1,Ipts)*CoefMu(1,Ipts) )*
     3                    BendFuncC(Ipts)
       ENDDO 


      ENDIF   !   What are ISurf and theta (Delta K) )


C
C       NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO RHO 
C
C       Value of the <A(B)|  DMCoeff(rho) |C>*Y_1^i * Y_3^j     -  DM = Dipole Moment 
C
         DipMomentBend(theta+1,ISurf,V2_1+1,V2_2+1,K1+1,NDeltaK,I0) =          
     1                 SimpsonIntegral(NumPtsRho,RhoMaximum,Phi_V_Phi)
        

      ENDDO    !   -   end of loop ever V2_2 quantum number 
      ENDDO    !   -   end of loop ever V2_1 quantum number 
      ENDIF    !   -   K2 <=KMAXC 
      ENDDO    !   -   end of loop ever DeltaK = K2-K1  
      ENDDO    !   -   end of loop ever K1  quantum number 
      ENDDO    !   -   theta 
      ENDDO    !   -   A,B - electronic states  Left  < |
      ENDDO    !   -   end of loop OVER ALL DIFERENT Y_1^i * Y_3^j COMBINATIONS


      CLOSE ( NFileMu )

C
      RETURN 

1009  WRITE (6,2009) NUNIT
2009  FORMAT(1H0,'ABTC_TOC.ERR  SCRATCH FILE COULD NOT BE O',
     1          'PENED, UNIT NUMBER IS ',I2)


      END 

C########################      BendMatrixElement   #########################
C###########################################################################


C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : READ FROM FILE COS(GAMMA_k) AND SIN(GAMMA_K) 
C             : RHO DEPENDENT FUNCTIONS 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE ReadCosSinGamma(CosGamma,SinGamma)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)
C
      CHARACTER*200 LINE
      CHARACTER*80 FILNAM
      INTEGER NP,K

      include 'dipolsys.h'
      include 'mu.h'

      REAL*8  CosGamma(KMAXAB+1,NumPtsRho),
     1        SinGamma(KMAXAB+1,NumPtsRho)



C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      WRITE(6,3001)
3001  FORMAT(1H0,'                 READ COS GAMMA AND SIN GAMMA '/)

     
C
C     OPEN FILE COS(GAMMA) AND SIN(GAMMA) FUNCTION FOR THE LOWER STATES
C
      
      FILNAM=DirInput(1:NumDirLetters)//FileGamma
      OPEN(UNIT=91,FILE=FILNAM)


C     SKIP FIRST LINES
      READ(91,1000) LINE
      READ(91,1000) LINE
1000  FORMAT(A20)

C
C     READ COS(GAMMA_k)- RHO DEPENDED FUNCTION - CosGamma
C

c      DO K =0,KMAXAB
        READ(91,*) ((CosGamma(K+1,NP), NP=1,NumPtsRho),K =0,KMAXAB)
c      ENDDO

C
C     READ COS(GAMMA_k)- RHO DEPENDED FUNCTION - CosGamma
C

        READ(91,*) ((SinGamma(K+1,NP), NP=1,NumPtsRho),K =0,KMAXAB)
c      DO K =0,KMAXAB
c       write(6,*) k,NumPtsRho
c        READ(91,*) (SinGamma(K+1,NP), NP=1,NumPtsRho)
c      ENDDO


      CLOSE(91)

      RETURN

      END
C#########################      ReadCosSinGamma     ########################
C###########################################################################



C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : SETS UP THE RHO DEPENDENT WAVE FUNCTIONS NEEDED TO
C             : CALCULATE THE MATRIX ELEMENTS OF THE DIPOLE MOMENT OF THE 
C             : TRIATOMIC MOLECULE USING THE MORBID SCHEME.
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE ReadBendFunction(NFileA,NFileB,NFileC)
C
      CHARACTER*200 LINE
      CHARACTER*80 FILNAM
      INTEGER NP,K,V2

      include 'dipolsys.h'
      include 'mu.h'

      REAL*8  BendFunc(NumPtsRho)

      INTEGER NFileA,NFileB,NFileC,NUNIT,KREC


C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      WRITE(6,3001)
3001  FORMAT(1H0,'                 READ BENDING WAVE FUNCTIONS.'/)


C
C     OPEN FILE WITH LOWER STATE BENDING FUNCTION
C

      FILNAM=DirInput(1:NumDirLetters)//FileBendAB
      OPEN(UNIT=92,FILE=FILNAM)

C

C
C     READ A- STATE BENDING BASIS FUNCTION - BendFuncA
C

C
C     SKIP FIRST LINES
C
      READ(92,1000) LINE
      READ(92,1000) LINE

1000  FORMAT(A20)

      DO K =0,KMAXAB
      DO V2=0,NV2Max(1)
c        DO NP=1,NumPtsRho
c           READ(92,*) BendFuncAB(1,V2+1,K+1,NP)
c        ENDDO
         KREC=K*( NV2Max(1)+1 )+V2+1
         READ(92,*) (BendFunc(NP), NP=1,NumPtsRho)
         WRITE (NFileA,REC=KREC)    (BendFunc(NP),NP=1,NumPtsRho)
      ENDDO
      ENDDO

C
C     READ B- STATE BENDING BASIS FUNCTION - BendFuncB
C

C
C     SKIP 2 FIRST LINES
C
      READ(92,1000) LINE
      READ(92,1000) LINE


      DO K =0,KMAXAB
      DO V2=0,NV2Max(2)
        READ(92,*) (BendFunc(NP), NP=1,NumPtsRho)

        KREC=K*( NV2Max(2)+1 )+V2+1
        WRITE (NFileB,REC=KREC)    (BendFunc(NP),NP=1,NumPtsRho)

      ENDDO
      ENDDO

      CLOSE(92)



C
C     OPEN FILE WITH UPPER STATE BENDING FUNCTION
C

      FILNAM=DirInput(1:NumDirLetters)//FileBendC
      OPEN(UNIT=92,FILE=FILNAM)
C
C     SKIP 2 FIRST LINES
C
      READ(92,1000) LINE
      READ(92,1000) LINE

C
C     READ C- STATE BENDING BASIS FUNCTION - BendFuncC
C

      DO K =0,KMAXC
      DO V2=0,NV2MaxC
        READ(92,*) (BendFunc(NP), NP=1,NumPtsRho)
        KREC=K*( NV2MaxC+1 )+V2+1
        WRITE (NFileC,REC=KREC)    (BendFunc(NP),NP=1,NumPtsRho)
      ENDDO
      ENDDO


      CLOSE(92)


      RETURN


      END

C#########################      ReadBendFunction    ########################
C###########################################################################


C##########################################################################
C#             :
C# DATE        : 05.06.2000
C# AUTHOR      : Yurchenko
C# UPDATES     :
C# LANGUAGE    : FORTRAN
C# PURPOSE     : NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT 
C# SUBPROGRAMS : 
C#      CALLED :
C#             :
C###########################################################################

      DOUBLE PRECISION FUNCTION  SimpsonIntegral(Npoints,XMax,F)  
C     IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I, Npoints

      REAL*8 Feven,Fodd,F0,Fmax,XMax,F(0:Npoints),h
C

C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      h = Xmax/Npoints            !   integration width  
      Feven=0                     
      Fodd =0
      F0   =F(0)
      Fmax =F(Npoints)

C
C     SUMMATION OF ODD AND EVEN PARTS OF THE SIMPSON INTEGRAL 
C

      DO I = 1,NPoints-2,2
           Fodd   = Fodd  + F(i  )
           Feven  = Feven + F(i+1)
      ENDDO 

      Fodd   = Fodd  + F(NPoints-1)

      SimpsonIntegral =  h/3.0D+00*
     1                 ( 4.0D+00*Fodd + 2.0D+00*Feven + F0 + Fmax)

      RETURN 
      END FUNCTION 

C#########################      SimpsonIntegral    #########################
C###########################################################################