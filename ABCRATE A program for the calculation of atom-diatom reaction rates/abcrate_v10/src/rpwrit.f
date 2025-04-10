!!!*************************************************************
! 文件/File: rpwrit.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpwrit.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPWRIT
C
C     RPWRIT - write out reaction path info
C
C  Called by:
C     RPHSUM - summarize reaction path info
C
C  Calls:
C     RPW2   - write out reaction path info
C     TITLES - print a 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *   NSSP(NSADDM), LRSTRT
      COMMON /RP2/    SLM, SLP
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *   FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *   X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      OPEN (UNIT=1, FILE='abc.1', FORM='FORMATTED',                     GCL1096
     *      STATUS='UNKNOWN', ERR=10)                                   GCL1096
      REWIND 1
      CALL TITLES (1, 1, 1)
      WRITE (1, 100) DELSV, DEL, ACASY, EPSASY, SLM, SLP, DELCUR,
     * DELCSP, SM, SP, THETA
      WRITE (1, 101) NSMAX, NSSP, NMID, IBOPT
      CALL RPW2 (NSMAX, SS, VS, DS, XKS, FBS, QFBS, GBS, XMOMS, X2,
     *   Y2, UXS, UYS, CAPS)
C
      RETURN
C
10    WRITE (6, 900)                                                    GCL1096
      STOP 'RPWRIT 1'                                                   GCL1096
  100 FORMAT (' DELSV,DEL,ACASY,EPSASY,SLM,SLP,DELCUR,DELCSP,',
     *   'SM,SP,THETA=', /, (1X, 1PE19.10, 3E20.10))
  101 FORMAT (' NSMAX,NSSP,NMID,IBOPT=', /, (1X, I9, 7I10))
900   FORMAT(/,1X,T5,'Error: Could not open abc.1 in RPWRIT')           GCL1096
      END
