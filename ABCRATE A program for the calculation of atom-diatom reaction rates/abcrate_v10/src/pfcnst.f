!!!*************************************************************
! 文件/File: pfcnst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pfcnst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PFCNST (BET, IS, S, D, XK, QVIB, CVIB, LP, LWKB)
C
C     PFCNST - compute stretching part. fcn.
C
C  Called by:
C     DATAIN - read in data and set up constants
C     PFCN   - compute GTS partition functions and free energy
C
C  Calls:
C     QMPART - compute part. fcn. for stretch
C     WKBINT - interpolate WKB energy levels from grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LP, LWKB
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
C
      IF (LSTATE .EQ. 0) THEN
         T = 2.0D0*D*BET/XK
         QVIBH = -T - LOG(1.D0-EXP(-2.D0*T))
         IF (LGS(3)) THEN
            CALL QMPART (BET, IS, S, D, XK, QVIB, LP, LWKB)
         ELSE
            QVIB = QVIBH
         END IF
      ELSE
         E = 2.D0*TNP1*D/XK
         QVIBH = -BET*E
         IF (LGS(3)) THEN
            E = E*(1.D0-0.5D0*TNP1/XK)
            IF (LWKB) CALL WKBINT (IS, S, NSTATE, E, UL, UG, PER)
            QVIB = -BET*E
         ELSE
            QVIB = QVIBH
         END IF
      END IF
      T = QVIB-QVIBH
      IF (T .GT. 80.D0) THEN                                            GCL1092
         WRITE (6, 6000) IS, LSTATE, LGS(3), LWKB, BET, S, D, XK, QVIB,
     *      QVIBH
         IS = 99999
      ELSE
         CVIB = EXP(T)
      END IF
      RETURN
 6000 FORMAT(2X,T5,'Warning: In PFNCST CVIB is greater than 1.E35',
     *     /,2X,T14,'IS, LSTATE,LGS3,LWKB=',2I5, 2L1,
     *     /,2X,T14,'BET,S,D,XK,QVIB,QVIBH=', 1P,6E13.5) 
      END
