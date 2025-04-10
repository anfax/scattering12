!!!*************************************************************
! 文件/File: rpread.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpread.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPREAD (ICHECK)
C
C     RPREAD - read in reaction path info
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     RPR2   - read in reaction path info
C     TITLES - print a 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TH(2)
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
      WRITE (6, 600)
C
      OPEN (UNIT=1, FILE='abc.1', FORM='FORMATTED',                     GCL1096
     *      STATUS='OLD', ERR=10)                                       GCL1096
      CALL TITLES (-1, 1, 2)
      CALL TITLES (1, 6, 2)
      READ (1, 100) DELSV, DEL, ACASY, EPSASY, SLM,
     * SLP, DELCUR, DELCSP, SM, SP, TH
      READ (1, 101) NSMAX, NSSP, NMID, IBOP
      CALL RPR2(NSMAX, SS, VS, DS, XKS, FBS, QFBS, GBS, XMOMS, X2,
     * Y2, UXS, UYS, CAPS)
      IF (ICHECK .NE. 0) THEN
         IF (IBOPT .NE. IBOP .OR. ABS(TH(1)-THETA(1)) .GT. 1.D-6 .OR.
     *      ABS(TH(2)-THETA(2)) .GT. 1.D-6) THEN
            WRITE(6, 6001) IBOP, IBOPT, TH, THETA
            STOP 'RPREAD 1'
         END IF
      ELSE
         THETA(1) = TH(1)
         THETA(2) = TH(2)
         IBOPT = IBOP
      END IF
      RETURN
10    WRITE (6, 900)                                                    GCL1096
      STOP 'RPREAD 2'                                                   GCL1096
900   FORMAT(/,1X,T5,'Error: Cannot open abc.1 in RPREAD')              GCL1096
  100 FORMAT (/, (1X, E19.10, 3E20.10))
  101 FORMAT (/, (1X, I9, 7I10))
  600 FORMAT (/,1X,T5,'The reaction coordinate information is read ',
     *                'in from the file',
     *        /,1X,T5,'that is linked to FORTRAN unit 1',
     *        /,1X,T5,'The title card from creation run is:')
 6001 FORMAT (/,1X,T5,'Error: The bend fit used in the saved data ',
     *                'is different from',
     *       /,1X,T12,'that specified by the new data.',
     *       /,1X,T12,'Old IBOPT=',I5,' new IBOPT=', I5,
     *       /,1X,T12,'Old THETA values=', F10.5,1X,F10.5,
     *       /,1X,T12,'New THETA values=', F10.5,1X,F10.5)
      END
