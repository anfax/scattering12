!!!*************************************************************
! 文件/File: pfcnb.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pfcnb.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PFCNB (BET, FB, AB, GB, QB, CB, DELPHI, LP)
C
C     PFCNB  - compute bending part. fcn.
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     DATAIN - read in data and set up constants
C     PFCN   - compute GTS partition functions and free energy
C
C  Calls:
C     COBPRT - compute partition function of centrifugal osc.
C     PRTQRT - compute partition function of harm-quartic potential
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LP
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/  LGS2
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
C
C Compute the harmonic partition function if FB > 0
      QBH = 0.D0                                                        GCL1096
      CB = 0.D0                                                         GCL1096
      IF (FB .GT. 0.D0) THEN                                            GCL1096
          WB = SQRT(FB*GB)                                              GCL1096
          IF (LSBEND .EQ. 0) THEN                                       GCL1096
              T = 0.5D0*WB*BET                                          GCL1096
              QBH = -T - LOG(1.D0-EXP(-2.D0*T))                         GCL1096
              QBH = QBH + QBH                                           GCL1096
          ELSE                                                          GCL1096
              EB = DBLE(NBND1+NBND2+1)*WB                               GCL1096
              QBH = -BET*EB + LOG(DGBND)                                GCL1096
          END IF                                                        GCL1096
          CB = 1.D35                                                    GCL1096
      ENDIF                                                             GCL1096
C
      IF (LGS2(6)) THEN
C  centrifugal oscillator partition function, delphi is not set 
         CALL COBPRT(BET, FB, AB, GB, QB, DELPHI)                       GCL0896
      ELSE IF (LGS(4)) THEN
C  uncoupled bending partition function
         CALL PRTQRT (BET, FB, AB, GB, QB, DELPHI, LP)
      ELSE
         DELPHI = 1000.D0                                               GCL1092
         QB = 80.D0                                                     GCL1092
         CB = 1.0D0
         IF (FB .GT. 0.D0) THEN
            QB = QBH                                                    GCL1096
            EVAL = MAX(1.D0/BET, 0.5D0*WB)                              GCL1096
            DELPHI = 2.D0*SQRT(2.D0*EVAL*GB/(WB*WB))                    GCL1096
         END IF
      END IF
C
      IF (FB .GT. 0 .AND. (LGS2(6) .OR. LGS(4))) THEN                   GCL1096
          DIFF = QB - QBH                                               GCL1096
          IF (DIFF .LT. 80.D0) CB = EXP(DIFF)                           TCA0497
      ENDIF                                                             GCL1096
C
      RETURN
      END
