!!!*************************************************************
! 文件/File: pmuomt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pmuomt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

         SUBROUTINE PMUOMT (ENGY, PROB, NQKPDM, NPETYP)
C
C        PSAG   - compute muOMT-type probabilites
C
C        Called by:
C                  KAPVA  - compute kappas
C        Calls:
C              nothing
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER (NQGKDM=81)
         COMMON /ESPEC/  ESPEC(40), NESPEC
         COMMON /PSAGCM/ NPSAG
         COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
         COMMON /QUADTH/ PT2(NQGKDM), WT2(NQGKDM,2), NQ22, NSEG2
         DIMENSION ENGY(NQKPDM), PROB(2, NQKPDM, NPETYP)
C
C   set the pointers for the muOMT, SCSAG, and LCG3 probabilities in the
C   array PROB (called PE in KAPVA).  
         PARAMETER (NMUOMT = 11)
         PARAMETER (NSCSAG =  3)
         PARAMETER (NLCG3  = 10)
         PARAMETER (IUNIT  = 6)
C
C   determine the number of energies at which the tunneling probabilities 
C   were calculated
         IF (NESPEC .LE. 0) THEN
             NE = NQ12
         ELSE
             NE = NESPEC
         END IF
         NE = NE*NSEG
C
C   initialize the variable NQ2
         NQ2 = (NQ22-1)/2
C
C   Loop over the two quadrature schemes.
C   The SCSAG probabilities are stored in PE in the following location
C   PROB(I,NQKPDM,3), and the LCG3 probabilities are stored in 
C   PROB(I,NQKPDM,10), where I = 1, 2 (the two different quadrature schemes 
C   used to determine the probabilities).  
C   The muOMT probabilites will be stored in PROB(I,NQKPDM,11)
C
         DO 20 I1 = 1, 2
               DO 10 I2 = 1, NE
                     PROB(I1,I2,NMUOMT) = MAX(PROB(I1,I2,NSCSAG),
     1                                        PROB(I1,I2,NLCG3))
10             CONTINUE
20       CONTINUE
C
C   If NPSAG = 0, write out the probabilities to the file linked to
C   FORTRAN unit IUNIT.
               IF (NPSAG .EQ. 0) THEN
                   WRITE (IUNIT, 600) NSEG2, NQ2, NSEG2, NQ22
                   WRITE (IUNIT, 610) 
     1             (ENGY(I), (PROB(J,I,NMUOMT),J=1,2),I=1,NE)
               ENDIF
C
600   FORMAT(/,2X,T5,'muOMT tunneling probabilities',
     1       /,2X,T10,'Energy (kcal)',T30,'N = ',I2,'*',I2,
     2                                    T50,'N = ',I2,'*',I2)
610   FORMAT(2X,T10,1PE13.5,T30,1PE13.5,T50,1PE13.5)
C
          RETURN
          END
