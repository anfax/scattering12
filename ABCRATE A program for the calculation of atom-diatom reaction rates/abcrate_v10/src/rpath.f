!!!*************************************************************
! 文件/File: rpath.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpath.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPATH
C
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C     REFLEC - reflect MEP info across symmetry point
C     RPREAD - read in reaction path info
C     RPSEG  - compute reaction path for a segment
C     RPHSUM - summarize reaction path info
C     SHIFT  - shift r.p. info in grid
C     TITLES - print out 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LRFLC, LTERM
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     *OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /RP2/    SLM, SLP
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE ISAD                                                         TCA1097
      DATA ISAD /0/                                                     GCL0795
C
      NSMAX = NSDM
      WRITE (6, 600)
C  restart option
      IF (LRSTRT) THEN
         CALL RPREAD (1)
      ELSE
         CALL TITLES (1, 6, 1)
      END IF
C  print out information about step sizes, etc.
      T1=  0.2D0*DELSV                                                  GCL1092
      T2 = 3.D0*T1                                                      GCL1092
      WRITE (6, 601) DELSV, T1, T2, DEL                                 GCL0795
      WRITE (6, 602) ACASY, EPSASY, SLM, SLP                            GCL0795
      WRITE (6, 603) DELCUR, DELCSP, DELSV
      WRITE (6, 604) SM, SP
      WRITE (6, 606) IBOPT
      IF (IBOPT .EQ. 3)                                                 GCL1092
     *    WRITE (6, 608) (180.D0*(1.D0 - THETA(I)/PI), I=1,2)           GCL1092
C  restart option not set, have to do some work
      IF (.NOT.LRSTRT) THEN
C  initialize some variables
         WRITE (6, 610)
         LRFLC = .FALSE.
         LTERM = .FALSE.
         IS = 0
         N = NSADMX(1)
         NSSP(1) = 1
         DELFIX = DEL
         ERREPS = .01D0*DEL                                             GCL1092
         SS(1) = 0.D0                                                   GCL1092
         VTEST = (2.D0*VSP(N) - VDELTA)/DBLE(NSDM)                      GCL1092
         XKTEST = ABS(2.D0*XKTS(N) - XKAB - XKBC)/DBLE(NSDM)            GCL1092
C  loop over saddle points
10       CONTINUE
            ISAD = ISAD + 1
            S = 0.D0                                                    GCL1092
C           WRITE (91, 9100) ISAD, NSAD
C9100       FORMAT (' RPATH: ISAD, NSAD=',2I5)
            IF (NSAD .GT. 0) THEN
C  ISGN = -1, go from saddle point towards reactants
               ISGN = -1
               ISMN = IS + 1
               ISMX = NSMAX
               IS = ISMX
C  IEND = -1,  end in reactant region
C       =  0,  end in intermediate well
               IF (ISAD .EQ. 1) THEN
                  IEND = -1
               ELSE
                  IEND = 0
               END IF
C              WRITE (91, 9101) ISGN, IEND, ISMN, ISMX, IS
C9101          FORMAT (' ISGN,IEND,ISMN,ISMX,IS=', 5I5)
               CALL RPSEG (DELFIX,ERREPS,VTEST,XKTEST,ISWL,IEND,ISAD,   GCL0493
     *                     ISGN,IS,ISMN,ISMX,LTERM,LRFLC)               GCL0493
C              WRITE (91, 9102) ISMN, ISMX, IS, S
C9102          FORMAT (' ISMN,ISMX,IS,S=', 3I5, 1PE13.5)
C  shift stored values so IS=ISMN
               IS = IS + 1
               ISHIFT = ISMN - IS
               IF (IEND .EQ. 0) THEN
                  SDIF = SS(ISWL) - S
               ELSE
                  SDIF = 0.D0                                           GCL1092
               END IF
               IF (ISHIFT .NE. 0) CALL SHIFT (IS, ISHIFT, ISGN, SDIF,
     *            ISAD, NSSP)
               ISMX = NSMAX
               WRITE (6, 612)
            END IF
C  reset IS to saddle point value for next segment, reflection, or
C  termination
            IS = NSSP(ISAD)
C  check for symmetric system
            LRFLC = LSYM .AND. ABS(R1SP(ISAD)-R2SP(ISAD)) .LT. 0.001D0  GCL1092
            IF (.NOT.(LTERM .OR. LRFLC)) THEN
               S = SS(IS)
C  ISGN = 1, go from saddle point towards products
               ISGN = 1
               ISMN = IS
               ISMX = NSMAX
C  IEND =  0,  end in intermediate well
C       =  1,  end in product region
               IF (ISAD .LT. NSAD) THEN
                  IEND = 0
               ELSE
                  IEND = 1
               END IF
C              WRITE (91, 9101) ISGN, IEND, ISMN, ISMX, IS
               CALL RPSEG (DELFIX,ERREPS,VTEST,XKTEST,ISWL,IEND,ISAD,   GCL0493
     *                     ISGN,IS,ISMN,ISMX,LTERM,LRFLC)               GCL0493
C              WRITE (91, 9102) ISMN, ISMX, IS, S
               IS = IS - 1
            END IF
         IF (ISAD.LT.NSAD .AND. .NOT.LTERM .AND. .NOT.LRFLC) GO TO 10
         NSMAX = IS
C  for symmetric system, reflect across symmetric S value
         IF (LSYM) THEN
            IF (LRFLC) THEN
               CALL REFLEC (IS, ISAD)
            ELSE
C  if LSYM is set but didn't hit the symmetric S value, write error
C     message and exit
               WRITE (6, 6000)
               STOP 'RPATH 1'
            END IF
         END IF
      END IF
C  summary section
      CALL RPHSUM
C  only return in subroutine, otherwise terminated by call exit
      RETURN
C
600   FORMAT (/, 1X, 6('*'), ' Reaction coordinate calculation',
     *        ' (all quantities in atomic units) ', 6('*'))
601   FORMAT (/,1X,T5, 'Reaction path information is stored every', 
     *        1PE13.5,  ' bohr until',
     *        /,1X,T5, 'the asymptotic region and at additional ',
     *                 'points', E13.5, /,1X,T5,'and', E13.5, 
     *                 ' bohr on either side of the saddle point.',
     *       //,1X,T5, 'The reaction coordinate is found in (X,Y) ',
     *                 'space by following',
     *        /,1X,T5, 'the gradient with an increment DEL = ', E13.5)
602   FORMAT(/,1X,T5, 'The system is considered in the asymptotic ',
     *                'region if the cosine',
     *       /,1X,T5, 'of the angle between the gradient and the ',
     *                'asymptotic gradient is',
     *       /,1X,T5, 'greater than', 0PF15.10,
     *       //,1X,T5, 'The calculations extend into the reactant and ',
     *                 'product regions until',
     *       /,1X,T5, 'the potential, Morse parameter, and bend ',
     *                'frequency are',
     *       /,1X,T5, 'within ', 1PE13.5, ' relative error compared ',
     *                'to the asymptotic values',
     *       /,1X,T5, 'or within the limits on s of ', 1PE13.5, 
     *                ' and ', 1PE13.5)
603   FORMAT (//,1X,T5, 'Step size for the curvature calculation is ',
     * 1PE13.5, /,1X,T5,'except for ', E13.5, ' within ', E13.5, 
     *                  ' bohr of the saddle point.',
     *        /,1X,T5, 'The value at the saddle point is ',
     *                 'interpolated from',
     *        /,1X,T5, 'one value on either side of saddle point.')
  604 FORMAT (/,1X,T5,'The curvature is calculated from ', 0PF10.6, 
     *                ' to', 0PF10.6)
  606 FORMAT (/,1X,T5, 'Option for fitting the bending potential ',
     *                 '(IBOPT) = ', I5)
  608 FORMAT (1X,T5,'The angles used for the fit are ', F10.5, 
     *              ' and ', F10.5, ' degrees.')
  610 FORMAT (78X, 'Gradient vector'/ 6X, 's', 10X, 'RAB', 8X, 'RBC',   TCA1097
     *   9X, 'X', 10X, 'Y', 10X, 'V', 13X, 'DELX/DEL', 5X, 'DELY/DEL',  TCA1097
     *   4X, 'Curvature', 2X, 'LASY', 2X, 'NSTAB'/)                     TCA1097
  612 FORMAT (1X)
 6000 FORMAT (/,2X,T5,'Error: The symmetry option was set but the ',
     *                'reflection point was not hit.')
      END
