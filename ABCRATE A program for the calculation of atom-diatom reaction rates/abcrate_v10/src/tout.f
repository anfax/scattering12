!!!*************************************************************
! 文件/File: tout.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tout.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TOUT (T, R, N)
C
C     TOUT   - print out table
C
C  Called by:
C     TABL20 - print out table of ratios
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LFLG
      CHARACTER*1 OUTPUT(132)
      DIMENSION R(6)
      DO 10 I = 1,132
         OUTPUT(I) = ' '
   10 CONTINUE
      IT = INT(T)
      REM = T - IT
      IF (REM .GE. .5D0) IT = IT + 1
      IC = 5
      LFLG = .TRUE.
      ITENS = 100000
   20 CONTINUE
         ITENS = ITENS/10
         I = IT/ITENS
         IF (.NOT.LFLG .OR. I .NE. 0) THEN
            OUTPUT(IC) = CHAR(I+48)
            LFLG = .FALSE.
         END IF
         IT = IT - I*ITENS
         IC = IC + 1
      IF (ITENS .GT. 1) GO TO 20
      IC0 = IC + 2
      ICP = IC0 + 4
      DO 100 I = 1,N
         R1 = R(I)
         IF (R1 .NE. 0.0D0) THEN
            IF (R1 .LE. 9950.0D0) THEN
               IC = IC0
               TENS = 1000.0D0
   30          CONTINUE
               IF (R1 .GT. 0.995D0*TENS .OR. IC .GT. IC0+6) GO TO 40
                  TENS = 0.1D0*TENS
                  IC = IC + 1
                  IF (IC .EQ. ICP) IC = IC + 1
               GO TO 30
   40          CONTINUE
               IF (R1 .GT. 0.995D0*TENS) THEN
                  NDIGIT = MIN(3, IC0+8-IC)
                  ITEN = 10**NDIGIT
                  R1 = R1 * .1D0 * ITEN / TENS
                  I1 = INT(R1)
                  REM = R1 - I1
                  IF (REM .GE. .5D0) I1 = I1 + 1
                  ISIG = 0
                  IF (IC .GE. ICP) THEN
                     ICT = ICP - 1
   50                CONTINUE
                        OUTPUT(ICT) = CHAR(48)
                        ICT = ICT + 1
                     IF (ICT .LT. IC) GO TO 50
                  END IF
   60             CONTINUE
                     ISIG = ISIG + 1
                     ITEN = ITEN/10
                     I2 = I1/ITEN
                     I1 = I1 - I2*ITEN
                     OUTPUT(IC) = CHAR(I2+48)
                     IC = IC + 1
                     IF (IC .EQ. ICP) IC = IC + 1
                  IF (ITEN .GT. 1) GO TO 60
                  OUTPUT(ICP) = '.'
                  IF (IC .LE. ICP) THEN
                     OUTPUT(IC) = CHAR(48)
                     IC = IC + 1
                  END IF
               ELSE
                  IC = IC0 + 3
                  ITEN = 1 - INT(LOG10(R1))                             GCL0895
                  IF (ITEN .GT. 35) THEN
                     OUTPUT(ICP) = CHAR(48)
                  ELSE
                     R1 = R1 * 10.0D0**ITEN
                     IT = INT(R1)
                     REM = R1 - IT
                     IF (REM .GE. .5D0) IT = IT + 1
                     IF (IT .GE. 10) THEN
                        ITEN = ITEN - 1
                        IT = 1
                     END IF
                     OUTPUT(IC) = CHAR(IT+48)
                     OUTPUT(ICP) = '.'
                     OUTPUT(IC+2) = 'E'
                     OUTPUT(IC+3) = '-'
                     IF (ITEN .LE. 9) THEN
                        OUTPUT(IC+4) = CHAR(ITEN+48)
                     ELSE
                        II = ITEN/10
                        OUTPUT(IC+4) = CHAR(II+48)
                        II = ITEN - II*10
                        OUTPUT(IC+5) = CHAR(II+48)
                     END IF
                  END IF
               END IF
            ELSE
               IC = IC0 + 3
               ITEN = INT(LOG10(R1))                                    GCL0895
               IF (ITEN .GE. 35) THEN
                  ITEN = 99
                  R1 = 9.1D0
               ELSE
                  R1 = R1 / 10.0D0**ITEN
               END IF
               IT = INT(R1)
               REM = R1 - IT
               IF (REM .GE. .5D0) IT = IT + 1
               IF (IT .GE. 10) THEN
                  ITEN = ITEN + 1
                  IT = 1
               END IF
               OUTPUT(IC) = CHAR(IT+48)
               OUTPUT(ICP) = '.'
               OUTPUT(IC+2) = 'E'
               OUTPUT(IC+3) = '+'
               IF (ITEN .LE. 9) THEN
                  OUTPUT(IC+4) = CHAR(ITEN+48)
               ELSE
                  II = ITEN/10
                  OUTPUT(IC+4) = CHAR(II+48)
                  II = ITEN - II*10
                  OUTPUT(IC+5) = CHAR(II+48)
               END IF
            END IF
         END IF
         IC0 = IC0 + 10
         ICP = ICP + 10
  100 CONTINUE
      WRITE(20, 2000) (OUTPUT(IC),IC=1,IC0)
      RETURN
 2000 FORMAT(10X, 132A1)
      END
