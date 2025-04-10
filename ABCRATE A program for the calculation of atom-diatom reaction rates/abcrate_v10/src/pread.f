!!!*************************************************************
! 文件/File: pread.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pread.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PREAD (L1,L2,NE,NQKPDM,EZ,VM,IOPTAU,ESV,               TCA0997
     *  ALFSV,ALFSV2,SLG,SRG,THET1,THETA,NTPSV,PERCNT,CSA0L,CSA0R,
     *  NFINAL,IEND,IE0,IC3D)
C
C     PREAD  - read in info needed for LAG calculation
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     TITLES - print a 3-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL L1,L2,L1T,L2T
      LOGICAL LEND,LERR
      DIMENSION ALFSV2(2,NQKPDM),ALFSV(2,NQKPDM),CSA0L(2,NQKPDM,4),
     *  CSA0R(2,NQKPDM,4),ESV(NQKPDM),PERCNT(2,2,NQKPDM),
     *  SLG(NQKPDM),SRG(NQKPDM),THETA(2,NQKPDM,4),THET1(2,NQKPDM),
     *  NTPSV(NQKPDM)
C
C  IE0 - last energy point successfully read in
C  IEND < 0, error in read
C       = 0, no eof encountered
C       = 1, eof encountered
      IU = IC3D + 10
C
      IF (IU .EQ. 11) THEN                                              GCL1096
          OPEN (UNIT=IU, FILE='abc.11', FORM='FORMATTED',               GCL1096
     *          STATUS='OLD', ERR = 900)                                GCL1096
      ELSE                                                              GCL1096
          OPEN (UNIT=IU, FILE='abc.12', FORM='FORMATTED',               GCL1096
     *          STATUS='OLD', ERR = 900)                                GCL1096
      ENDIF                                                             GCL1096
      NFTSV = -1
      IESV = -1
      IE = 0
      REWIND IU
      CALL TITLES (-1,IU,3)
      LEND = .FALSE.
      LERR = .FALSE.
      READ (IU,1000,ERR=10,END=11) EZT,VMT,L1T,L2T,NET,IOP2             TCA0997
         IEND = 0
         GO TO 15
   10 CONTINUE
         LERR = .TRUE.
         GO TO 15
   11 CONTINUE
         LEND = .TRUE.
   15 CONTINUE
      IF (LERR .OR. LEND .OR. ABS(EZT-EZ) .GT. 1.D-4 .OR.
     *   ABS(VMT-VM) .GT. 1.D-4 .OR. NE .NE. NET .OR.                   TCA0997
     *   IOPTAU .NE. IOP2) THEN                                         TCA0997
         IEND = -1
         WRITE (6,6001) EZT,VMT,L1T,L2T,NET,IOP2,EZ,VM,                 TCA0997
     *      L1,L2,NE,IOPTAU
      ELSE
   20    CONTINUE
            IE = 0
            READ (IU,1001,ERR=25,END=26) NFT
               GO TO 27
   25       CONTINUE
               LERR = .TRUE.
               GO TO 50
   26       CONTINUE
               LEND = .TRUE.
               GO TO 50
   27       CONTINUE
            IF (NFT .GT. NFINAL) THEN
               IEND = -2
               WRITE (6,6002) NFT,NFINAL
               GO TO 50
            END IF
            IF (NFT .EQ. NFINAL) THEN
               WRITE (6,601) IU
               CALL TITLES (1,6,3)
            END IF
   30       CONTINUE
               IE = IE + 1
               READ(IU,1002,ERR=40,END=41) ESV(IE),NTPSV(IE),
     *          (ALFSV(I,IE),I=1,2),(ALFSV2(I,IE),I=1,2),
     *          SLG(IE),SRG(IE),(THET1(I,IE),I=1,2),
     *          ((THETA(I,IE,J),I=1,2),J=1,4),
     *          ((PERCNT(I,J,IE),I=1,2),J=1,2),
     *          ((CSA0L(I,IE,J),I=1,2),J=1,4),
     *          ((CSA0R(I,IE,J),I=1,2),J=1,4)
                  GO TO 42
   40          CONTINUE
                  LERR = .TRUE.
                  GO TO 42
   41          CONTINUE
                  LEND = .TRUE.
   42          CONTINUE
            IF (.NOT.LERR .AND. .NOT.LEND .AND. IE .LT. NE) GO TO 30
            IF (LERR .OR. LEND) IE = IE - 1
            NFTSV = NFT
            IESV = IE
         IF (NFT. NE. NFINAL) GO TO 20
   50    CONTINUE
      END IF
      IE0 = IE
      IF ((LERR .OR. LEND) .AND. IEND .EQ. 0) THEN
         IF (NFT .EQ. NFINAL .OR. (NFTSV .EQ. NFINAL-1 .AND. IESV .EQ.
     *      NE)) THEN
            IEND = 1
         ELSE
            IEND = -3
            WRITE (6,6003) NFT,NFINAL
         END IF
      END IF
      RETURN
C
900   WRITE (6, 901) IU                                                 GCL1096
      STOP 'PREAD 1'                                                    GCL1096
901   FORMAT(/,1X,T5,'Error: Cannot open the LAG file linked to ',      GCL1096
     *               'FORTRAN unit ', I5)                               GCL1096
600   FORMAT (/,2X,T5,'Probability data read in from file',2X,A8,'*',A8,
     *        /,2X,T5,'Title card from creation run:')
601   FORMAT (/,2X,T5,'Probability data read in from the file ',
     *                'attached to unit',I3,
     *        /,2X,T5,'Title card from creation run:')
1000  FORMAT(/,1X,1P,2E20.10,2(2X,L1),3I10) 
1001  FORMAT(/,1X,I10)
1002  FORMAT(1X,1PE20.10,I10,/,(1X,1P,6E20.10))  
6001  FORMAT (/,2X,T5,'Problem with LAG probability file, old and new',
     * ' data do not agree',/,10X,'EZ',11X,'VM',12X,'LGS(9)',5X,
     * 'LGS(10)',10X,'NEMAX',5X,'IOPTAU',/,' OLD',                      TCA0997
     * 1P,2E13.5,10X,L1,10X,L1,10X,I5,8X,I5,6X,I5,/,' NEW', 
     * 2E13.5,10X,L1,10X,L1,10X,I5,8X,I5)                               TCA0997
6002  FORMAT (/,2X,T5,'Problem with LAG probability file, the smallest',
     * ' value of NFINAL from the restart file is',I5,
     * /,2X,T5,'which is larger than the value requested',I5)
6003  FORMAT (/,2X,T5,'Problem with LAG probability file, the last ',
     *  'value of NFINAL from the restart file is',I5,
     *  /,2X,T5,'which is smaller than the value requested',I5)
      END
