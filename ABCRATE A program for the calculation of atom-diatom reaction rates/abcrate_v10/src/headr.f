!!!*************************************************************
! 文件/File: headr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: headr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

C
C***********************************************************************
C  HEADR
C***********************************************************************
C
      SUBROUTINE headr
C
C   This subprogram prints out the ABCRATE program header to unit 6.
C
C   Called by:
C             VTSTM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      WRITE (6,1000)
      WRITE (6,1100)
      WRITE (6,1200)
      WRITE (6,1300)
      WRITE (6,1400)
      WRITE (6,1500)
      WRITE (6,1000)
C
C
      RETURN
C
 1000 FORMAT (/,1X,T2,78('*'))
 1100 FORMAT (/,2X,T24,'ABCRATE-version 10.0 (October 1997)',
     *        //,2X,T40,'by')              
 1200 FORMAT (/,2X,T33,'Bruce C. Garrett',
     *        /,2X,T24,'Molecular Science Research Center',
     *        /,2X,T23,'Pacific Northwest National Laboratory',
     *        /,2X,T31,'Richland, Washington')
 1300 FORMAT (/,2X,T12,'Gillian C. Lynch, Thomas C. Allison, and ',
     *                 'Donald G. Truhlar',
     *        /,2X,T16,'Department of Chemistry and Supercomputer ',
     *                 'Institute',/,2X,T29,'University of Minnesota',
     *        /,2X,T30,'Minneapolis, Minnesota')
 1400 FORMAT (///,2X,T15,'A program for the calculation of rate ',
     *                   'constants by', /,2X,T15,'generalized ',
     *                   'transition state theory (GTST) with ',
     *          /,2X,T15,'semiclassical tunneling probabilities for',
     *          /,2X,T15,'collinear-dominated atom-diatom reactions.')
 1500 FORMAT (///,2X,T15,'Copyright, 1996, Donald G. Truhlar and ',
     *          /,2X,T15,'Regents of the University of Minnesota')
C
       END
