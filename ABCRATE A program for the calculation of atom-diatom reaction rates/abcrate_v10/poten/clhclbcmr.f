!!!*************************************************************
! 文件/File: clhclbcmr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: clhclbcmr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

C
      SUBROUTINE PREPOT
C
C   System:          ClHCl
C   Functional form: LEPS (London-Eyring-Polanyi-Sato)
C   Common name:     BCMR
C   Reference:       D. K. Bondi, J. N. L. Connor, J. Manz, and J. Romelt
C                    J. Mol. Phys. 50, 467 (1983)
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohrs
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state 
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Cl
C      infinitely far from the HCl diatomic and R(HCl) set equal to the
C      HCl equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(first Cl-H)
C                            R(2) = R(H-second Cl)
C                            R(3) = R(first Cl-second Cl)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, energy, and derivatives of
C        the energy with respect to the coordinates.
C      /PT2CM/ NDER, NDUM(21)
C        passes the control flags where
C        NDER  = 0 => no derivatives are calculated
C        NDER  = 1 => calculate first derivatives of the energy for the 
C                     ground electronic state with respect to the coordinates
C        NDUM  - not used 
C      /PT4CM/ IPTPRT 
C        passes the FORTRAN unit number used for POTENTIAL output
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A+BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy
C        of the C atom "infinitely" far from the AB diatomic and R(AB) set 
C        equal to Re(AB), the equilibrium bond length for the AB diatomic.  
C        In this potential the AB valley represents Cl infinitely far from
C        the ClH diatomic and R(ClH) equal to Re(ClH).
C        Similarly, the terms EASYBC and EASYAC represent the energies in the 
C        other HCl and the ClCl asymptotic valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NDER             1
C      IPTPRT           6
C
C*****
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER (R2     = 1.41421356D0)
         PARAMETER (CKCAU  = 627.5095D0)
         PARAMETER (CANGAU = 0.52917706D0)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
C
      DIMENSION ZPO(3), OP3Z(3), ZP3(3),TZP3(3), TOP3Z(3), DO4Z(3), B(3)
      DIMENSION X(3),COUL(3),EXCH(3)
C
         WRITE (IPTPRT, 600) D, RE, BETA, Z
C
C   CONVERT TO ATOMIC UNITS
      DO  10 I = 1,3
      D(I)=D(I)/CKCAU 
      RE(I) = RE(I)/CANGAU
      BETA(I) = BETA(I)*CANGAU
C   COMPUTE USEFUL CONSTANTS
      ZPO(I) = 1.0D0 + Z(I)
      OP3Z(I) = 1.0D0 + 3.0D0*Z(I)
      TOP3Z(I) = 2.0D0*OP3Z(I)
      ZP3(I) = Z(I) + 3.0D0
      TZP3(I) = 2.0D0*ZP3(I)
      DO4Z(I) = D(I)/4.0D0/ZPO(I)
   10 B(I) = BETA(I)*DO4Z(I)*2.0D0
C   Initialize the energy in the three asymptotic valleys
      EASYAB = D(1)
      EASYBC = D(2)
      EASYAC = D(3)
C
      RETURN
C
      ENTRY POT
C
C   System:          ABC
C   Functional form: LEPS (London-Erying-Polyani-Sato)
C   Reference:       S. Sato
C                    J. Chem. Phys. 592, 2465 (1955)
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
         ENDIF
C
      DO 20 I = 1,3
            X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
            COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
            EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
20    CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+
     1      (EXCH(3)-EXCH(1))**2)
      ENERGY = COUL(1) + COUL(2) + COUL(3) - (RAD/R2) + D(2)
C
C   Compute the derivatives of the LEPS energy w.r.t.the coordinates
      IF (NDER .EQ. 1) THEN
          SVAL = EXCH(1) + EXCH(2) + EXCH(3)
          DO 30 I = 1,3
                DEDR(I)=0.0D0 
                IF(X(I).LT.1.0D-30) GO TO 30
                TVAL = (3.0D0*EXCH(I)-SVAL)/R2*(OP3Z(I)*X(I)-ZP3(I)) 
                IF(ABS(RAD) .LT. 1.0D-32 .AND. 
     1             ABS(TVAL) .GT. 1.0D-12) THEN 
                   WRITE(IPTPRT, 6000) TVAL, RAD
                ELSE IF(ABS(RAD).GT.1.0D-32) THEN 
                   TVAL = TVAL/RAD 
                END IF 
                DEDR(I)=B(I)*X(I)*(TVAL-ZP3(I)*X(I)+OP3Z(I))  
30        CONTINUE
      ENDIF
C
600   FORMAT(/,1X,'*****',1X,'Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'ClHCl BCMR LEPS potential energy surface',
     *       //,1X,T5,'Parameters:',
     *        /,1X,T5,'Bond', T46, 'Cl-H', T58, 'H-Cl', T68, 'Cl-Cl',
     *        /,1X,T5,'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /,1X,T5,'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
6000  FORMAT(/,1X,'In the ClHCl potential BCMR, T, RAD = ',1P,2E15.7,
     1       /,1X,'T/RAD has been set equal to T')  
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT
         COMMON /PT2CM/ NDER, NDUM(21)
         COMMON /SATOCM/ D(3), RE(3), BETA(3), Z(3) 
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER /6, 1/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA D    / 106.477D0, 106.477D0, 57.983D0/
         DATA RE   / 1.275D0, 1.275D0, 1.988D0/    
         DATA BETA / 1.868D0, 1.868D0, 2.002D0/
         DATA Z    / 0.115D0, 0.115D0, 0.115D0/
C
         END
C*****
