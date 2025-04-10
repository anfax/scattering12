!!!*************************************************************
! 文件/File: rpstor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpstor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPSTOR(ISAD,ISGN,IS,ISMN,ISMX,IEND,LTERM)
C
C     RPSTOR - store rp info and check grid
C
C  Called by:
C     RPSEG  - compute reaction path for a segment
C
C  Calls:
C     SHIFT  - shift r.p. info in grid
C     STORE  - store reaction path info
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LTERM
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
C
C     WRITE (91, 9110) IS, S, ISMN, ISMX, NSMAX, NSSP(ISAD)
C9110 FORMAT (' SAVE IS,S,ISMN,ISMX,NSMAX,NSSP=', I5, F10.5, 4I5)
C  store rp info
      CALL STORE (IS)
C  increment grid location IS
      IS = IS + ISGN
      IF (IS .GT. ISMX .OR. IS .LT. ISMN) THEN
C  IS outside allowed range, try to shift or terminate
         IS = IS - ISGN
C           WRITE (91, 9111) IS, ISGN
C9111       FORMAT (' IS OUT OF RANGE, IS,ISGN=', 2I5)
         IF (ISGN .LT. 0 .AND. IEND .LE. 0) THEN
C  overran grid going towards reactants, try reducing ISMN before
C     shifting
            ISMN = ISMN - 1
            IF (ISMN .LE. 0) THEN
               WRITE (6, 6007)
               STOP 'RPSTOR 1'
            END IF
         END IF
         IF (ISGN .EQ. 1 .AND. IEND .EQ. 1 .AND. NSSP(ISAD) .LT.
     *    NSMAX/2) THEN
C              WRITE (91, 9112) NSSP(ISAD)
C9112          FORMAT (' SADDLE POINT CENTERED, NSSP=', I5)
            LTERM = .TRUE.
         ELSE
            IF (IEND .EQ. 1 .AND. ISGN .LT. 0) THEN
C  last segment, from asymptote towards well, shift IS-NSMAX one to the
C     right
               ISHIFT = 1
               II = -1
            ELSE
C  all other segments, shift 1-IS one to the left
               ISHIFT = -1
               II = 1
            END IF
            SDIF = 0.D0                                                 GCL1092
            CALL SHIFT (IS, ISHIFT, II, SDIF, ISAD, NSSP)
C              WRITE (91, 9113) ISHIFT, II, NSSP
C9113          FORMAT (' ISHIFT,II,NSSP=', 6I5)
            IF (NSSP(1) .LT. 1) THEN
               WRITE (6, 6008)
               STOP 'RPSTOR 2'
            END IF
         END IF
      END IF
      RETURN
 6007 FORMAT (/,1X,T5,'Error: In RPSTOR ran out of space in the grid, ',
     *                'ISMN=0')
 6008 FORMAT (/,1X,T5,'Error: in RPSTOR ran out of space in the grid, ',
     *                'first saddle point at IS=1.')
      END
