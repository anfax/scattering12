!!!*************************************************************
! 文件/File: kaptyp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: kaptyp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE KAPTYP
C
C     KAPTYP - set up labels and indices for kappas
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /KINDX/  NKP, NKT, IOB, NLAG
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32,NSEG3, IOPTAU,   TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
C
      DO 5 I=1,NKALAB                                                   TCA0197
    5   LPTYP(I)=.TRUE.                                                 TCA0197

      NKP = 4
      NLAG = NKP + 1
      KLAG = 0
      IF (LLAG) KLAG = 8                                                GCL93
      NKP = NKP + KLAG
      NKT = NKP + 3
      IOB = NKT + 3
C
      DO 10 I = 1,NKT
         ATYP(I) = 'GARBAG'
   10 CONTINUE
C
      ATYP(1) = 'TST/CAG'                                               TCA0197
      ATYP(2) = 'CVT/CAG'
      ATYP(3) = 'MEPSAG'
      ATYP(4) = 'MCPTSAG'
      LPTYP(4) = .FALSE.                                                TCA0197
      ATYP(5) = 'SCTSAG'
      LPTYP(5) = .FALSE.                                                TCA0197
      ATYP(6) = 'PATSAG'
      LPTYP(6) = .FALSE.                                                TCA0197
C
      ATYP(NKT) = 'MCPSAG'
      LPTYP(NKT) = .FALSE.                                              TCA0197
      ATYP(NKT+1) = 'CD-SCSAG'
      ATYP(NKT+2) = 'PASAG'
      LPTYP(NKT+2) = .FALSE.                                            TCA0197
C
      IF (LLAG) THEN
        ATYP(NLAG+2) = 'RLCG'
        LPTYP(NLAG+2) = .FALSE.                                         TCA0197
        ATYP(NLAG+3) = 'RLAG'
        LPTYP(NLAG+3) = .FALSE.                                         TCA0197
        ATYP(NLAG+4) = 'LAG'
        ATYP(NLAG+5) = 'LAG 0.0'
        LPTYP(NLAG+5) = .FALSE.                                         TCA0197
        ATYP(NLAG+6) = 'LAG 0.5'
        LPTYP(NLAG+6) = .FALSE.                                         TCA0197
        ATYP(NLAG+7) = 'LCG3'
        ATYP(NLAG+8) = 'muOMT'                                          GCL93
        ATYP(NLAG+9) = 'COMT'                                           GCL93
        LPTYP(NLAG+9) = .FALSE.                                         TCA0997
      END IF
      RETURN
      END
