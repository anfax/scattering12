!!!*************************************************************
! 文件/File: chkrst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: chkrst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C 
      SUBROUTINE CHKRST
C
C  CHKRST - checks that restart files exist
C
C  Called by:
C     DATAIN - read in data and set up constants
C
C  Calls:
C     nothing
C
      DOUBLE PRECISION DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY,
     *                 PT3, WT3 
      INTEGER NSMAX, NSSP, NQ32, NSEG3, IOPTAU, NQGKDM, NSADDM
      LOGICAL LRSTRT, LLAG, LLAGRS, LGS(10), LGS2(10), LEXIT
      PARAMETER (NQGKDM=81)
      PARAMETER (NSADDM=4)
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     * NSSP(NSADDM), LRSTRT
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,
     * LLAG, LLAGRS 
      COMMON /LOGIC/ LGS
      COMMON /LOGIC2/ LGS2
      SAVE LEXIT                                                        TCA1097
      DATA LEXIT /.FALSE./
C
      IF (LRSTRT) THEN
          OPEN (UNIT=1, FILE='abc.1', FORM='FORMATTED', STATUS='OLD',
     *          ERR=10)
          CLOSE (UNIT=1, STATUS='KEEP')
      ENDIF
          
11    IF (LGS2(5)) THEN
          OPEN (UNIT=14, FILE='abc.14', FORM='FORMATTED', STATUS='OLD',
     *          ERR=20)
          CLOSE (UNIT=14, STATUS='KEEP')
      ENDIF
12    IF (LLAGRS .AND. LGS(9)) THEN
          OPEN (UNIT=11, FILE='abc.11', FORM='FORMATTED', 
     *          STATUS='OLD', ERR=30)
          CLOSE (UNIT=11, STATUS='KEEP')
      ENDIF
13    IF (LLAGRS .AND. LGS(10)) THEN
          OPEN (UNIT=12, FILE='abc.12', FORM='FORMATTED', 
     *          STATUS='OLD', ERR=40)
          CLOSE (UNIT=12, STATUS='KEEP')
      ENDIF
C
      GO TO 15 
C
10    LEXIT = .TRUE.
      WRITE (6, 900)
      GO TO 11
20    LEXIT = .TRUE.
      WRITE (6, 910)
      GO TO 12
30    LEXIT = .TRUE.
      WRITE (6, 920)
      GO TO 13
40    LEXIT = .TRUE.
      WRITE (6, 930)
C
15    IF (LEXIT) STOP 'CHKRST 1'
C
900   FORMAT(/,1X,T5,'Error: Cannot open abc.1 file')
910   FORMAT(/,1X,T5,'Error: Cannot open abc.14 file')
920   FORMAT(/,1X,T5,'Error: Cannot open abc.11 file')
930   FORMAT(/,1X,T5,'Error: Cannot open abc.12 file')
C
      RETURN
      END
C
