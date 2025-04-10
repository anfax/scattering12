!!!*************************************************************
! 文件/File: pwrite.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pwrite.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PWRITE (L1,L2,NE,NDIM,EZ,VM,IOPTTT,                    TCA0997
     *ESV,ALFSV,ALFSV2,SLG,SRG,THET1,THETA,NTPSV,PERCNT,
     *CSA0L,CSA0R,NFINAL,IE,IC3D)
C
C     PWRITE - write out info needed to restart LAG calculation
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     TITLES - print a 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL L1,L2,LFRST1,LFRST2
      CHARACTER*3 STAT
      DIMENSION ALFSV(2,NDIM),ALFSV2(2,NDIM),CSA0L(2,NDIM,4),
     *CSA0R(2,NDIM,4),ESV(NDIM),PERCNT(2,2,NDIM),SLG(NDIM),
     *SRG(NDIM),THETA(2,NDIM,4),THET1(2,NDIM)
      DIMENSION NTPSV(NDIM)
      PARAMETER (NQGKDM=81,NQKPDM=4*NQGKDM)
      LOGICAL LLAG,LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM),WT3(NQGKDM,2), NQ32, NSEG3,IOPTAU,
     * LLAG,LLAGRS
      SAVE LFRST1, LFRST2                                               TCA1097
      DATA LFRST1/.TRUE./, LFRST2/.TRUE./
C
      IF (IC3D .EQ. 1) THEN
         IU = 11
         IF (LFRST1) THEN
            STAT = 'NEW'
            IF (LLAGRS) STAT = 'OLD'
            LFRST1 = .FALSE.
         END IF
      ELSE
         IU = 12
         IF (LFRST2) THEN
            STAT = 'NEW'
            IF (LLAGRS) STAT = 'OLD'
            LFRST2 = .FALSE.
         END IF
      END IF
      IF (STAT .EQ. 'NEW') THEN
          IF (IU .EQ. 11) THEN                                          GCL1096
              OPEN (UNIT=IU, FILE='abc.11', FORM='FORMATTED',           GCL1096
     *              STATUS='NEW', ERR=10)                               GCL1096
         ELSE                                                           GCL1096
              OPEN (UNIT=IU, FILE='abc.12', FORM='FORMATTED',           GCL1096
     *              STATUS='NEW', ERR=10)                               GCL1096
         ENDIF                                                          GCL1096
         CALL TITLES (1,IU,1)
         WRITE (IU,1000) EZ,VM,L1,L2,NE,IOPTTT                          TCA0997
      END IF
      IF (IE .EQ. 1) WRITE (IU,1001) NFINAL
      WRITE (IU,1002) ESV(IE),NTPSV(IE),
     *   (ALFSV(I,IE),I=1,2),(ALFSV2(I,IE),I=1,2),SLG(IE),SRG(IE),
     *   (THET1(I,IE),I=1,2),((THETA(I,IE,J),I=1,2),J=1,4),
     *   ((PERCNT(I,J,IE),I=1,2),J=1,2),
     *   ((CSA0L(I,IE,J),I=1,2),J=1,4),((CSA0R(I,IE,J),I=1,2),J=1,4)
      STAT = 'OLD'
      RETURN 
10    WRITE (6, 900) IU                                                 GCL1096
      STOP 'PWRITE 1'                                                   GCL1096
900   FORMAT(/,1X,T5,'Error: Cannot open the LAG file linked to ',      GCL1096
     *               'FORTRAN unit ', I5)                               GCL1096
1000  FORMAT(' EZ,VM,L1,L2,NE,IOPTTT=',/,1X,1P,2E20.10,2(2X,L1),2I10)   TCA0997
1001  FORMAT(' NFT=',/,1X,I10)
1002  FORMAT(1X,1PE20.10,I10,/,(1X,1P,6E20.10))                         GCL0992
      END
