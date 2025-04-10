!!!*************************************************************
! 文件/File: mutens.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: mutens.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE MUTENS ( DMMACY , DMMACZ , DMMBCX )
      IMPLICIT REAL*8 (A-H,O-Z)
C             :
C DATE        : 28.04.1986
C AUTHOR      : PER JENSEN
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATION OF THE EXPANSION COEFFICINTETS (IN THE
C             : S I'S, THE DELTA R I S, AND THE Y I'S) OF THE NON-
C             : VANISHING MU TENSOR ELEMENTS.
C SUBPROGRAMS :
C      CALLED : AMAT, TRNSSR, INVERT, TRNSRY
C             :
C INPUT       : RHO  INSTANTANEOUS VALUE OF THE BENDING COORDINATE
C             : EPS  INSTANTANEOUS VALUE OF THE ANGLE EPSILON DEFINED
C             :      IN P. JENSEN, COMP. PHYS. REP. 1, 1-55 (1983)
C             :
C OUTPUT      : MUABY(1:15) (AB = XX, YY, ZZ, RR, YZ, AND RX). EXPANSION
C             :             COEFFICIENTS OF THE MUAB TENSOR ELEMENT IN
C             :             TERMS OF THE YI PARAMETERS.
C             :
C             :
C    THE ARRAYS CONTAIN THE COEFFICIENTS OF THE FOLLOWING
C    QUANTITIES:
C
C    I             IABS(I)        IABR(I)           MUABY(I)
C  -------------------------------------------------------------
C    1             1              1                 1
C    2             S1             DR1               Y1
C    3             S3             DR3               Y3
C    4             S1**2          DR1**2            Y1**2
C    5             S3**2          DR3**2            Y3**2
C    6             S1*S3          DR1*DR3           Y1*Y3
C    7                            DR1**3            Y1**3
C    8                            DR3**3            Y3**3
C    9                            DR1**2*DR3        Y1**2*Y3
C   10                            DR1*DR3**2        Y1*Y3**2
C   11                            DR1**4            Y1**4
C   12                            DR3**4            Y3**4
C   13                            DR1**3*DR3        Y1**3*Y3
C   14                            DR1*DR3**3        Y1*Y3**3
C   15                            DR1**2*DR3**2     Y1**2*Y3**2
C
C
C
      INTEGER V1 ,V2, V3 , V2P1,
     1       NTEST 
C
      REAL*8 DNYS(6),DNZS(6),CMYS(6),CMZS(6),V3YS(6),V3ZS(6),
     1       DNYR(15),DNZR(15),DNYY(15),DNZY(15),U3YR(15),U3ZR(15),
     2       REC32(15),V3ZR(15),V3YR(15),V1YS(6),V1ZS(6),
     2       REC12(15),V1ZR(15),V1YR(15),U1YR(15),U1ZR(15),
     3       CRRC(15),AUXIY(15),DMY(15),DMZ(15),
     4       AUXIZ(15),UPYR(15),UPZR(15),
     3       DMPARR(15),DMPERP(15),DMYR(15),DMZR(15),DMMPXR(15)
C


      REAL*8  DIPBCX(53)

c ---- Quantum number parameters for stretching wave functions ----- 
      include 'dipolsys.h'

      COMMON /DIPBC/ DIPBCX 


      REAL*8 IXXS(6),IYYS(6),IZZS(6),IRRS(6),IYZS(6),IXRS(6),
     1       IXXR(15),IYYR(15),IZZR(15),IRRR(15),IYZR(15),IXRR(15),
     2       MUXXR(15),MUYYR(15),MUZZR(15),MURRR(15),
     3       MUYZR(15),MUXRR(15),
     4       MUXXY(15),MUYYY(15),MUZZY(15),MURRY(15),
     3       MUYZY(15),MUXRY(15),
     5       PX1R(15),PX3R(15),PR1R(15),PR3R(15),PX1Y(15),
     6       PX3Y(15),PR1Y(15),PR3Y(15),
     7       PXMXX1(15),PXMXX3(15),PRMXR1(15),PRMXR3(15),
     8       PXMXR1(15),PXMXR3(15),PRMRR1(15),PRMRR3(15),
     6       MUZZIN(15),MUYZIN(15)
      REAL*8  A1Y1,A1Z1,A2Y1,A2Z1,A3Y1,A3Z1,
     2       A1Y3,A1Z3,A2Y3,A2Z3,A3Y3,A3Z3,
     5       YE1,ZE1,YE2,ZE2,YE3,ZE3
      REAL*8  AMMA11,AMMA33,AMMA13,AMMA31,
     1       APMA11,APMA33,APMA13,APMA31,G11,G33,G13,G31,
     2       CHX11,CHX13,CHX31,CHX33,CHR11,CHR13,CHR31,CHR33
C
      include 'dipole.h'
C
      REAL*8 DMMACY(15) , DMMACZ(15) , DMMBCX(15) 
C
      REAL*8   AUXILI(15),
     1       SBARS(6),SBARR(15),AUXIR(15),AUXI2(15)
C
      INTEGER ZC1,ZC2,ZC3,ZELEC,II
      REAL*8 EC

      REAL*8 TMREF

      COMMON /CHARGE/ EC,ZC1,ZC2,ZC3,ZELEC
C
C

      include 'value.h'
      include 'molcul.h'
c      include 'dimen.h'
      include 'bcoeff.h'


C
C UPON INPUT OF RHO AND EPS, THE SUBROUTINE AMAT CALCULATES THE
C NECESSARY VALUES OF THE A MATRIX ELEMENTS, THE A PRIME
C MATRIX ELEMENTS, THE A VECTOR COMPONENTS AND THEIR FIRST
C AND SECOND DERIVATIVES.
C
      CALL AMAT ( A1Y1,A1Z1,A2Y1, A2Z1, A3Y1, A3Z1,
     2           A1Y3,A1Z3,A2Y3,  A2Z3, A3Y3, A3Z3,
     5           YE1,ZE1,YE2,ZE2,YE3,ZE3)
C
C
C     CALCULATE CENTER-OF-MASS COORDINATES FOR THE REFERENCE
C     MOLECULE (EXPANSION COEFFIENTS IN THE S EXPANSION)
C
      TMREF=REFM1+REFM2+REFM3
      CMYS(1)=(REFM1*YE1 + REFM2*YE2 + REFM3*YE3)/TMREF
      CMYS(2)=(REFM1*A1Y1 + REFM2*A2Y1 + REFM3*A3Y1)/TMREF
      CMYS(3)=(REFM1*A1Y3 + REFM2*A2Y3 + REFM3*A3Y3)/TMREF
      CMYS(4)=0.0D+00
      CMYS(5)=0.0D+00
      CMYS(6)=0.0D+00
      CMZS(1)=(REFM1*ZE1 + REFM2*ZE2 + REFM3*ZE3)/TMREF
      CMZS(2)=(REFM1*A1Z1 + REFM2*A2Z1 + REFM3*A3Z1)/TMREF
      CMZS(3)=(REFM1*A1Z3 + REFM2*A2Z3 + REFM3*A3Z3)/TMREF
      CMZS(4)=0.0D+00
      CMZS(5)=0.0D+00
      CMZS(6)=0.0D+00
C
C     CALCULATE NUCLEAR DIPOLE MOMENT COORDINATES
C     (EXPANSION COEFFIENTS IN THE S EXPANSION)
C
C
      DNYS(1)=(ZC1*YE1 + ZC2*YE2 + ZC3*YE3
     .        -ZELEC*CMYS(1))*EC
      DNYS(2)=(ZC1*A1Y1 + ZC2*A2Y1 + ZC3*A3Y1
     .        -ZELEC*CMYS(2))*EC
      DNYS(3)=(ZC1*A1Y3 + ZC2*A2Y3 + ZC3*A3Y3
     .        -ZELEC*CMYS(3))*EC
      DNYS(4)=0.0D+00
      DNYS(5)=0.0D+00
      DNYS(6)=0.0D+00
      DNZS(1)=(ZC1*ZE1 + ZC2*ZE2 + ZC3*ZE3
     .        -ZELEC*CMZS(1))*EC
      DNZS(2)=(ZC1*A1Z1 + ZC2*A2Z1 + ZC3*A3Z1
     .        -ZELEC*CMZS(2))*EC
      DNZS(3)=(ZC1*A1Z3 + ZC2*A2Z3 + ZC3*A3Z3
     .        -ZELEC*CMZS(3))*EC
      DNZS(4)=0.0D+00
      DNZS(5)=0.0D+00
      DNZS(6)=0.0D+00
C
C     TRANSFORM TO POLYNOMIALS IN THE DELTA R'S
C
      CALL TRNSSR ( DNYS , DNYR )
      CALL TRNSSR ( DNZS , DNZR )
C
C     TRANSFORM TO POLYNOMIALS IN THE Y'S
C
      CALL TRNSRY ( DNYR , DNYY )
      CALL TRNSRY ( DNZR , DNZY )
C
C     CALCULATE COORDINATES FOR 2-1 VECTOR
C     (EXPANSION COEFFIENTS IN THE S EXPANSION)
C
C
      V1YS(1)=YE1-YE2
      V1YS(2)=A1Y1-A2Y1
      V1YS(3)=A1Y3-A2Y3
      V1YS(4)=0.0D+00
      V1YS(5)=0.0D+00
      V1YS(6)=0.0D+00
      V1ZS(1)=ZE1-ZE2
      V1ZS(2)=A1Z1-A2Z1
      V1ZS(3)=A1Z3-A2Z3
      V1ZS(4)=0.0D+00
      V1ZS(5)=0.0D+00
      V1ZS(6)=0.0D+00
C
C     TRANSFORM TO POLYNOMIALS IN THE DELTA R'S
C
      CALL TRNSSR ( V1YS , V1YR )
      CALL TRNSSR ( V1ZS , V1ZR )
C
C     PREPARE EXPANSION COEFFICIENTS FOR 1/R12
C
      REC12( 1)=1.0D+00/RE12
      REC12( 2)=-1.0D+00/RE12**2
      REC12( 3)=0.0D+00
      REC12( 4)=1.0D+00/RE12**3
      REC12( 5)=0.0D+00
      REC12( 6)=0.0D+00
      REC12( 7)=-1.0D+00/RE12**4
      REC12( 8)=0.0D+00
      REC12( 9)=0.0D+00
      REC12(10)=0.0D+00
      REC12(11)=1.0D+00/RE12**5
      REC12(12)=0.0D+00
      REC12(13)=0.0D+00
      REC12(14)=0.0D+00
      REC12(15)=0.0D+00
C
C     CALCULATE COORDINATES OF UNIT VECTOR IN THE
C     2-1 DIRECTION. THE EXPANSION COEFFICIENTS ARE
C     CONTAINED IN THE ARRAYS U1YR AND U1ZR
C
      CALL PLMLT2 ( REC12 , V1YR , U1YR )
      CALL PLMLT2 ( REC12 , V1ZR , U1ZR )
C
C     CALCULATE COORDINATES FOR 2-3 VECTOR
C     (EXPANSION COEFFIENTS IN THE S EXPANSION)
C
C
      V3YS(1)=YE3-YE2
      V3YS(2)=A3Y1-A2Y1
      V3YS(3)=A3Y3-A2Y3
      V3YS(4)=0.0D+00
      V3YS(5)=0.0D+00
      V3YS(6)=0.0D+00
      V3ZS(1)=ZE3-ZE2
      V3ZS(2)=A3Z1-A2Z1
      V3ZS(3)=A3Z3-A2Z3
      V3ZS(4)=0.0D+00
      V3ZS(5)=0.0D+00
      V3ZS(6)=0.0D+00
C
C     TRANSFORM TO POLYNOMIALS IN THE DELTA R'S
C
      CALL TRNSSR ( V3YS , V3YR )
      CALL TRNSSR ( V3ZS , V3ZR )
C
C     PREPARE EXPANSION COEFFICIENTS FOR 1/R32
C
      REC32( 1)=1.0D+00/RE32
      REC32( 2)=0.0D+00
      REC32( 3)=-1.0D+00/RE32**2
      REC32( 4)=0.0D+00
      REC32( 5)=1.0D+00/RE32**3
      REC32( 6)=0.0D+00
      REC32( 7)=0.0D+00
      REC32( 8)=-1.0D+00/RE32**4
      REC32( 9)=0.0D+00
      REC32(10)=0.0D+00
      REC32(11)=0.0D+00
      REC32(12)=1.0D+00/RE32**5
      REC32(13)=0.0D+00
      REC32(14)=0.0D+00
      REC32(15)=0.0D+00
C
C     CALCULATE COORDINATES OF UNIT VECTOR IN THE
C     2-3 DIRECTION. THE EXPANSION COEFFICIENTS ARE
C     CONTAINED IN THE ARRAYS U3YR AND U3ZR
C
      CALL PLMLT2 ( REC32 , V3YR , U3YR )
      CALL PLMLT2 ( REC32 , V3ZR , U3ZR )
C
C     IF (REFM1 .EQ. REFM3) THEN
      IF (DABS(REFM1-REFM3).LT.0.000001) THEN
C
C     GENERATE EXPANSION COEFFIENTS FOR 1/SQRT(2*(COS(RHOBAR)+1)
C     (CONTAINED IN CRRC)
C
      CALL RECCOS ( CRRC )
C
C     GENERATE COORDINATES FOR A UNIT VECTOR IN THE 'PARALLEL'
C     DIRECTION
C
      DO 5000 II=1,15
5000  AUXIY(II) = U3YR(II)-U1YR(II)
      CALL PLMLT2 ( CRRC , AUXIY , UPYR )
C
      DO 5100 II=1,15
5100  AUXIZ(II) = U3ZR(II)-U1ZR(II)
      CALL PLMLT2 ( CRRC , AUXIZ , UPZR )
C
      ELSE
C
      DO 5105 II=1,15
      UPYR(II) = U3YR(II)
5105  UPZR(II) = U3ZR(II)
C
      ENDIF
C
C
C
C
C     GENERATE EXPANSION COEFFIENTS FOR THE ELECTRONIC
C     DIPOLE MOMENT FUNCTIONS (PARALLEL AND PERPENDICULAR)
C     TAKEN FROM INPUT.
C
      CALL GENDIP ( DMPARR , DMPERP , DMMPXR )
C
C     GENERATE EXPANSION COEFFIENTS FOR THE ELECTRONIC
C     DIPOLE MOMENT ALONG THE Y AXIS.
      CALL PLMLT2 ( DMPARR , UPYR , AUXIY )
      CALL PLMLT2 ( DMPERP , UPZR , AUXIZ )
      DO 5200 II=1,15
5200  DMYR(II) = AUXIY(II)+AUXIZ(II)
C
C     GENERATE EXPANSION COEFFIENTS FOR THE ELECTRONIC
C     DIPOLE MOMENT ALONG THE Z AXIS.
C
      CALL PLMLT2 ( DMPARR , UPZR , AUXIY )
      CALL PLMLT2 ( DMPERP , UPYR , AUXIZ )
      DO 5300 II=1,15
5300  DMZR(II) = AUXIY(II)-AUXIZ(II)
C
C     TRANSFORM TO POLYNOMIALS IN THE Y'S
C
      CALL TRNSRY ( DMYR , DMY )
      CALL TRNSRY ( DMZR , DMZ )
      CALL TRNSRY (DMMPXR,DMMBCX)
C
C     ADD NUCLEAR PART
C
      DO II=1,15
         DMMACY(II) = DMY(II)+DNYY(II)
         DMMACZ(II) = DMZ(II)+DNZY(II)
      ENDDO  


C      DO 5400 II=1,15
C        IF (NUMBER.EQ.1) THEN
C          DMMMY(II) = DMY(II)+DNYY(II)
C          DMMMZ(II) = DMZ(II)+DNZY(II)
C        ELSE
C          DMPPY(II) = DMY(II)+DNYY(II)
C          DMPPZ(II) = DMZ(II)+DNZY(II) 
C        ENDIF  
C5400  CONTINUE


      RETURN
      END
C
C
