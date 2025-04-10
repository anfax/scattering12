!!!*************************************************************
! 文件/File: vspln2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vspln2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VSPLN2 (IC3D,NFINAL,NSMAX,SMAX,VPN)
C
C     VSPLN2 - spline fit of adiabatic potential in product channel
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     EBEND  - compute bending energy levels
C     LOCS  - locate position of S in grid such that SS(IS) < or = S <
C        SS(IS+1)
C     SPL1D1 - spline fit
C     SPL1B1 - spline fit
C     WKB    - compute WKB energy levels for stretch
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IOP(2)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/  LGS2
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /OPTION/ IOPT(20)
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      COMMON /SPLNV2/ VV2(NSDM1), AV2(NSDM1), BV2(NSDM1), CV2(NSDM1),
     *DV2(NSDM1), SCRTCH(NSDM1), NS2, ISN
      COMMON /TP2CM/ ULS(NSDM1),UGS(NSDM1)
      COMMON /VIBTA2/ TAUS(NSDM1)
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      NS2 = NSMAX
C Get location of smax in grid
      CALL LOCS(ISN,SMAX)
C Compute excited state adiabatic potential along reaction coordinate
      TNP1 = 2.0D0*NFINAL + 1.D0                                        GCL1092
      DO 30 IS = 1, NS2
         S = SS(IS)
         D = DS(IS)
         XK = XKS(IS)
         E = 2.D0 * D * TNP1 / XK
         IF (LGS(3) .OR. LGS(7)) THEN
            E = E * (1.D0 - 0.5D0 * TNP1/XK)
            IF(LGS(7)) CALL WKB(IS,S,NFINAL,E,ULS(IS),UGS(IS),TAUS(IS))
         END IF
         VV2(IS) = VS(IS) + E
         IF (IC3D .EQ. 2) THEN
            IF(LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
               CALL COBEND(NBND1,NBND2,FBS(IS),QFBS(IS),GBS(IS),EB)
            ELSE
C  uncoupled bending energy level
               EB = EBEND(NBND1, FBS(IS), QFBS(IS), GBS(IS))
               IF (NBND1 .EQ. NBND2) THEN
                  EB = 2.D0*EB
               ELSE
                  EB = EB + EBN(NBND2, FBS(IS), QFBS(IS), GBS(IS))
               END IF
            END IF
            VV2(IS) = VV2(IS) + EB
         END IF
   30 CONTINUE
C Spline fit potential curve
      IOP(1) = 5
      IOP(2) = 5
      CALL SPL1D1(NS2,SS,VV2,SCRTCH,IOP,1,AV2,BV2,CV2)
      CALL SPL1B1(NS2,SS,VV2,SCRTCH,1,AV2,BV2,CV2,DV2)
C     WRITE(6, 6603) NN, AV2(1), BV2(1), CV2(1), DV2(1)
C6603 FORMAT(' NN=', I5, /' COEFF', 2X, 1P4E15.7)
      IF (IOPT(5) .EQ. 0) THEN
         E = 2.D0 * DAB * TNP1 / XKAB
         IF (LGS(3) .OR. LGS(7)) THEN
            E = E * ( 1.D0 - 0.5D0 * TNP1 / XKAB)
            T = 50.0D0
            IF (LGS(7)) CALL WKB (-2, T, NFINAL, E, UL, UG, PER)
         END IF
         VPN = VDELTA + E
      ELSE
         VPN = VV2(NS2)
      END IF
C     WRITE (6, 6620) (TAUS(IS),IS=1,NS2)
C6620 FORMAT (' TAUS='/ (1X, 1P8E15.7))
      RETURN
      END
