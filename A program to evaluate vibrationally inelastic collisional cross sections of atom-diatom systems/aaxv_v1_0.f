!!!*************************************************************
! 文件/File: aaxv_v1_0.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: aaxv_v1_0.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

! AAXVEDWAVE.  A PROGRAM TO EVALUATE VIBRATIONALLY INELASTIC COLLISIONAL  AAXV0000
! 1   CROSS SECTIONS OF ATOM-DIATOM SYSTEMS.  M.M. NOVAK.                 AAXV0000
! REF. IN COMP. PHYS. COMMUN. 46 (1987) 417                               AAXV0000
               PROGRAM EDWAVE                                           AAXV0001
C**************************************************************         AAXV0002
C**************************************************************         AAXV0003
C                                                                       AAXV0004
C            V(E)DW-CALCULATIONS ON ATOM-DIATOM SYSTEMS                 AAXV0005
C                   M.M. NOVAK (1984)                                   AAXV0006
C            IMPLEMENTED IN FORTRAN 77 M.M. NOVAK (1986)                AAXV0007
C                 based on a program by                                 AAXV0008
C  G.G.BALINT-KURTI, J.H.V.LENTHE, L.ENO, R.SAKTREGER (1979)            AAXV0009
C**************************************************************         AAXV0010
C**************************************************************         AAXV0011
C                                                                       AAXV0012
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0013
C                                                                       AAXV0014
      COMMON CR(8000)                                                   AAXV0015
      COMMON/CORE / KSCR,KSR,KSI,KETA,KCRJ,KTOTAL,KXK,KFK,KCRJG,        AAXV0016
     #                        NCHM,NCORE,KCR,KRR,KPOTRX,KDPOTRX         AAXV0017
C                                                                       AAXV0018
C NCORE SHOULD BE EQUAL TO THE DIMENSION OF ARRAY CR IN                 AAXV0019
C        BLANK COMMON (SEE SUBR. 'MAINCS')                              AAXV0020
C                                                                       AAXV0021
      NCORE = 8000                                                      AAXV0022
C                                                                       AAXV0023
      CALL INITIALIZE                                                   AAXV0024
C                                                                       AAXV0025
      CALL EDW(CR(KCR+KSCR),CR(KCR+KSR),CR(KCR+KSI),CR(KCR+KETA),       AAXV0026
     #         CR(KCR+KETA+NCHM),CR(KCR+KETA+NCHM*2),                   AAXV0027
     #         CR(KCR+KETA+NCHM*3),CR(KCR+KXK),CR(KCR+KFK),             AAXV0028
     #         CR(KCR+KETA+NCHM*4),CR(KCR+KETA+NCHM*5),                 AAXV0029
     #         CR(KCR+KCRJG),CR(KCR+KCRJ),CR(KCR+KTOTAL),               AAXV0030
     &         CR(KCR+KRR),CR(KCR+KPOTRX),CR(KCR+KDPOTRX))              AAXV0031
C                                                                       AAXV0032
      END                                                               AAXV0033
C                                                                       AAXV0034
      BLOCK DATA                                                        AAXV0035
C---------------------------------------------------------------------- AAXV0036
C DATA INITIALIZATION FOR EDWAVE                                        AAXV0037
C                                                                       AAXV0038
C MEANING OF THE COMMON BLOCK VARIABLES :                               AAXV0039
C                                                                       AAXV0040
C /OPT/                                                                 AAXV0041
C     IPRNT  : PRINT-FLAG, DETERMINES THE QUANTITY OF OUTPUT            AAXV0042
C     IDWINTG: FLAG IF DISTORTED WAVE INTEGRALS SHOULD BE COMPUTED      AAXV0043
C     IDISC  : IF NE.0 ==> S-MATRICES ARE WRITTEN TO UNIT IDISC         AAXV0044
C                                                                       AAXV0045
C                                                                       AAXV0046
C /CHAN/                                                                AAXV0047
C     RMIN/RMAX : MINIMUM/MAXIMUM DISTANCE FOR DIST.WAVE INTEGRATION    AAXV0048
C     DFACTOR : USED IN FINDING START FOR INTEGRATION (SUBROUTINE START)AAXV0049
C     ACC   : ACCURACY-PARAMETER (COMPUTED FROM EPSIL IN SUBR. INITIALIZAAXV0050
C     NCHAN : NUMBER OF CHANNELS FOR A CROSS SECTION                    AAXV0051
C     MAXSTEP : MAXIMUM NUMBER OF STEPS ALLOWED FOR FINDING A           AAXV0052
C                                   TURNING POINT IN SUBROUTINE TURN    AAXV0053
C     KJ    : NUMBER OF OPEN CHANNELS  (CALCULATED IN SUBROUTINE ELEVEL)AAXV0054
C                                                                       AAXV0055
C /QUANT/                                                               AAXV0056
C     JTOT  : TOTAL ANGULAR MOMENTUM                                    AAXV0057
C     BRED  : 2*(REDUCED MASS)/HBAR**2                                  AAXV0058
C     ERED  : TOTAL ENERGY USED (SEE SUBR.EDW)                          AAXV0059
C     VIBQUAN: QUANTUM OF VIBRATIONAL ENERGY                            AAXV0060
C     IEDW  : SELECTS DISTORTED (-2) OR EXPONENTIAL DISTORTED (-1) WAVE AAXV0061
C     LWKBJ : DETERMINES IF WKBJ APPROXIMATION IS TO BE USED            AAXV0062
C     STWKBJ: CRITERIUM FOR START OF WKBJ APPROXIMATION                 AAXV0063
C                                                                       AAXV0064
C /PASS/                                                                AAXV0065
C     ETOT/ESTEP/NETOT : MINIMUM/STEPSIZE/NUMBER OF STEPS  FOR LOOP     AAXV0066
C                                                   OVER TOTAL ENERGY   AAXV0067
C     JJTMAX/JJTMIN/JTSTEP : MAXIMUM/MINIMUM/STEPSIZE FOR TOTAL         AAXV0068
C                                                    ANGULAR MOMENTUM   AAXV0069
C     CONV/CONVP : CONVERSION FACTORS FOR TOTAL AND PARTIAL J CROSSECTIOAAXV0070
C                                   RESPECTIVELY (SUBROUTINE WRITEC)    AAXV0071
C     NVMIN/NVMAX : MAXIMUM/MINIMUM FOR VIBRATIONAL QUANTUM NUMBER      AAXV0072
C                                                                       AAXV0073
C /CORE /                                                               AAXV0074
C     Kxxx  : STARTING ADDRESS OF ARRAY xxx IN BLANK COMMON             AAXV0075
C             COMPUTED IN SUBROUTINE INITIALIZE / USED IN MAIN PROGRAM  AAXV0076
C     NCHM  : MAXIMUM NUMBER OF CHANNELS SIMULTANEOUSLY NEEDED          AAXV0077
C     NCORE : NUMBER OF WORDS AVAILABLE IN BLANK COMMON                 AAXV0078
C                   (ONLY USED IN MAIN PROGRAM AND ROUTINE              AAXV0079
C                         MAINCS FOR SIMULATING DYNAMICAL CORE-USAGE)   AAXV0080
C     KCR   : EFFECTIVE STARTING POINT IN BLANK COMMON (ARRAY CR)       AAXV0081
C                                                                       AAXV0082
C /SURF /                                                               AAXV0083
C     IPNT  : POINTER USED WITH POTENTIAL MATRIX ELEMENTS               AAXV0084
C     NPOINT: NUMBER OF POINTS IN PROPAGATION GRID                      AAXV0085
C     NVMXP1: EQUAL TO NVMAX + 1                                        AAXV0086
C                                                                       AAXV0087
C /ANG  /                                                               AAXV0088
C     DRR   : PROPAGATION STEP SIZE                                     AAXV0089
C     NGAM  : NUMBER OF SCATTERING ANGLES                               AAXV0090
C     NROINT: NUMBER OF POINTS IN PROPAGATION GRID                      AAXV0091
C     IOLD  : POINTER - IS POTENTIAL MATRIX AVAILABLE?                  AAXV0092
C                                                                       AAXV0093
C /PARAM/                                                               AAXV0094
C     RVDW,DVDW,ALVDW: PARAMETERS FOR VAN DER WAALS MOLECULE            AAXV0095
C     RDIA,DDIA,ALDIA: PARAMETERS FOR DIATOMIC MOLECULE                 AAXV0096
C     ZXK   : PARAMETER DEFINED IN SUBROUTINE INITIALIZE                AAXV0097
C---------------------------------------------------------------------- AAXV0098
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0099
C                                                                       AAXV0100
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV0101
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV0102
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV0103
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV0104
     #              NVMIN,NVMAX                                         AAXV0105
      COMMON/CORE / KSCR,KSR,KSI,KETA,KCRJ,KTOTAL,KXK,KFK,KCRJG,        AAXV0106
     #              NCHM,NCORE,KCR,KRR,KPOTRX,KDPOTRX                   AAXV0107
C                                                                       AAXV0108
      DATA IPRNT, IDWINTG, IDISC /1,1,0/                                AAXV0109
      DATA DFACTOR, MAXSTEP /30.0D0,250/                                AAXV0110
      DATA IEDW /-1/                                                    AAXV0111
      DATA LWKBJ, STWKBJ /1,0.01D0/                                     AAXV0112
      DATA NETOT, JTSTEP, CONV, CONVP /1,10,-1.0D0,-1.0D0/              AAXV0113
      DATA NVMIN, NVMAX /0,1/                                           AAXV0114
      END                                                               AAXV0115
C                                                                       AAXV0116
      SUBROUTINE INITIALIZE                                             AAXV0117
C*******************************************************************    AAXV0118
C INPUT AND CORE-PARTITIONING FOR EDWAVE                                AAXV0119
C                                                                       AAXV0120
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0121
      COMMON CR(8000)                                                   AAXV0122
      DIMENSION TITLE(10)                                               AAXV0123
      DIMENSION D(8), E(8), W(8), V(8), Q(8), O(8), S(8)                AAXV0124
C                                                                       AAXV0125
      COMMON/CORE / KSCR,KSR,KSI,KETA,KCRJ,KTOTAL,KXK,KFK,KCRJG,        AAXV0126
     #              NCHM,NCORE,KCR,KRR,KPOTRX,KDPOTRX                   AAXV0127
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV0128
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV0129
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV0130
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV0131
     #              NVMIN,NVMAX                                         AAXV0132
      COMMON/SURF / IPNT,NPOINT,NVMXP1                                  AAXV0133
      COMMON/ANG  / DRR,NGAM,NROINT,IOLD                                AAXV0134
      COMMON/PARAM/ RVDW,DVDW,ALVDW,RDIA,DDIA,ALDIA,ZXK                 AAXV0135
C                                                                       AAXV0136
      NAMELIST/EDW/ IPRNT,IDWINTG,IDISC,MAXSTEP,JTMIN,JTMAX,JTSTEP,EMIN,AAXV0137
     #                ETOTAL,NESTEP,CONV,CONVP,HBAR,RSCALE,VSCALE,      AAXV0138
     #                RMASS,DFACTOR,RBEG,REND,EPSIL,IEDW,LWKBJ,STWKBJ,  AAXV0139
     $                VIBQUAN,NVMIN,NVMAX,NGAM,NROINT,                  AAXV0140
     $                RVDW,DVDW,ALVDW,RDIA,DDIA,ALDIA,DIMASS,IOLD       AAXV0141
C                                                                       AAXV0142
C INITIALIZE HEADING                                                    AAXV0143
C                                                                       AAXV0144
      DATA D/6HDDDD  ,6HDD DD ,4*6HDD  DD,6HDD DD ,6HDDDD  /            AAXV0145
      DATA E/6HEEEEEE,2*6HEE    ,2*6HEEEEEE,2*6HEE    ,6HEEEEEE/        AAXV0146
      DATA W/3*7HWW   WW,4*7HWW W WW,7H WWWWW /                         AAXV0147
      DATA V/5*6HVV  VV,2*6H VVVV ,6H  VV  /                            AAXV0148
      DATA Q/6HIIIIII,6*6H  II  ,6HIIIIII/                              AAXV0149
      DATA O/6H  OO  ,2*6H OOOO ,2*6HOO  OO,2*6H OOOO ,6H  OO  /        AAXV0150
      DATA S/6H  SS  ,6HSS  SS,6HSS    ,6H  SS  ,6H  SS  ,6H    SS,     AAXV0151
     *      6HSS  SS,6H  SS  /                                          AAXV0152
C                                                                       AAXV0153
C INITIALIZE CONSTANTS/VARIABLES                                        AAXV0154
C                                                                       AAXV0155
      DATA JTMIN, JTMAX /0,100/                                         AAXV0156
      DATA EMIN, ETOTAL, NESTEP /-1.0D0,-1.0D0,0/                       AAXV0157
      DATA HBAR, RSCALE, VSCALE /3*1.0D0/                               AAXV0158
      DATA RMASS /-1.0D0/, VIBQUAN /-1.0D0/                             AAXV0159
      DATA RBEG, REND /0.5D0,12.0D0/, EPSIL, AFACT /0.002,1.0D-4/       AAXV0160
      DATA MAXC /0/, NBYTE /8/                                          AAXV0161
      DATA NGAM, NROINT, IOLD /8,148,1/                                 AAXV0162
C                                                                       AAXV0163
C READ INPUT                                                            AAXV0164
C                                                                       AAXV0165
      OPEN (UNIT=5,FILE='DIA.DAT',STATUS='OLD')                         AAXV0166
      READ(5,500) TITLE                                                 AAXV0167
      READ(5,EDW)                                                       AAXV0168
C                                                                       AAXV0169
      ACC  = AFACT*EPSIL                                                AAXV0170
      BRED = 2.0D0*RMASS*VSCALE*(RSCALE/HBAR)**2                        AAXV0171
      RMIN = RBEG/RSCALE                                                AAXV0172
      RMAX = REND/RSCALE                                                AAXV0173
      JJTMAX = JTMAX + 1                                                AAXV0174
      JJTMIN = JTMIN + 1                                                AAXV0175
      IF (CONV.LT.0.0)    CONV = RSCALE**2                              AAXV0176
      IF (CONVP.LT.0.0)   CONVP = CONV                                  AAXV0177
      ESTEP   = 0.0D0                                                   AAXV0178
      ETOT    = ETOTAL/VSCALE                                           AAXV0179
      IF (EMIN.GE.0.0D0 .AND. NESTEP.NE.0) THEN                         AAXV0180
        NETOT = NESTEP + 1                                              AAXV0181
        ESTEP = (ETOTAL-EMIN)/(DFLOAT(NESTEP)*VSCALE)                   AAXV0182
        ETOT  = EMIN/VSCALE                                             AAXV0183
      END IF                                                            AAXV0184
      ZXK = DSQRT(2.0D0*DIMASS*DDIA)/(ALDIA*HBAR)                       AAXV0185
      DRR = (RMAX-RMIN+1.0D0)/NROINT                                    AAXV0186
C                                                                       AAXV0187
      K   = 100                                                         AAXV0188
      KK  = K*2 + 2                                                     AAXV0189
      KCR = MAINCS(CR,KK,NBYTE)                                         AAXV0190
      IF (KCR.LE.0) THEN                                                AAXV0191
        WRITE(6,601) KCR,KK                                             AAXV0192
        STOP 1                                                          AAXV0193
      END IF                                                            AAXV0194
C                                                                       AAXV0195
C RESERVE SPACE/DETERMINE CORE PARTITIONING                             AAXV0196
C                                                                       AAXV0197
      NVMXP1 = NVMAX + 1                                                AAXV0198
      NMAX   = NROINT*NVMXP1*(NVMXP1+1)/2                               AAXV0199
      NLEV   = (NVMAX-NVMIN) + 2                                        AAXV0200
      NCHM   = NLEV                                                     AAXV0201
C                                                                       AAXV0202
      NCHM   = 0                                                        AAXV0203
      JB     = NVMIN+1                                                  AAXV0204
      JE     = NVMAX+1                                                  AAXV0205
      DO 21 JJ = JB,JE                                                  AAXV0206
21    NCHM   = NCHM+JJ                                                  AAXV0207
      NCHM   = NCHM + 1                                                 AAXV0208
      NLSCR  = NCHM*(NCHM+1)/2                                          AAXV0209
c     IF (IABS(IEDW).EQ.2)      NLSCR = 0                               AAXV0210
C                                                                       AAXV0211
C COMMON BLOCK CORE                                                     AAXV0212
C                                                                       AAXV0213
      KTOTAL =  1                                                       AAXV0214
      KCRJ   = KTOTAL + NLEV*(NLEV+1)                                   AAXV0215
      KXK    = KCRJ + NLEV*(NLEV+1)/2                                   AAXV0216
      KFK    = KXK + NLEV                                               AAXV0217
      KSCR   = KFK + NLEV                                               AAXV0218
      KCRJG  = KSCR                                                     AAXV0219
      KSR    = KSCR + NLSCR                                             AAXV0220
      KSI    = KSR + NLSCR                                              AAXV0221
      KETA   = KSI + NLSCR                                              AAXV0222
      KRR    = KETA + 6*NCHM +15                                        AAXV0223
      KPOTRX = KRR + 150                                                AAXV0224
      KDPOTRX= KPOTRX + NMAX                                            AAXV0225
C                                                                       AAXV0226
      MAXC   = KDPOTRX + NMAX                                           AAXV0227
      MAXC   = MAX0(MAXC,KK)                                            AAXV0228
      KCR    = MAINCS(CR,MAXC,NBYTE)                                    AAXV0229
      IF (KCR.LE.0)      IPRNT = MAX0(IPRNT,10)                         AAXV0230
C                                                                       AAXV0231
C PRINT INPUT-PARAMETERS/CORE-PARTITIONING                              AAXV0232
C                                                                       AAXV0233
      IF (IEDW.NE.-1 .AND. IEDW.NE.-2)     IEDW = -1                    AAXV0234
      WRITE(6,605)                                                      AAXV0235
      DO 50 I = 1,8                                                     AAXV0236
      IF (IEDW.EQ.-1) WRITE(6,607) V(I),E(I),D(I),W(I),Q(I),O(I),S(I)   AAXV0237
      IF (IEDW.EQ.-2) WRITE(6,637) V(I),D(I),W(I),Q(I),O(I),S(I)        AAXV0238
   50 CONTINUE                                                          AAXV0239
      WRITE(6,608)                                                      AAXV0240
      WRITE(6,611) TITLE                                                AAXV0241
C                                                                       AAXV0242
      WRITE(6,612) NVMIN,NVMAX                                          AAXV0243
      IF (NESTEP.GT.0) WRITE(6,613) EMIN, ETOT, ETOTAL, NESTEP          AAXV0244
      IF (NESTEP.EQ.0) WRITE(6,614) ETOTAL, ETOT                        AAXV0245
      WRITE(6,615) JTMIN, JTMAX, JTSTEP                                 AAXV0246
      WRITE(6,617) HBAR, RSCALE, VSCALE                                 AAXV0247
      WRITE(6,624) CONV, CONVP                                          AAXV0248
      WRITE(6,618) RMASS, BRED                                          AAXV0249
      WRITE(6,619) RBEG, RMIN, REND, RMAX, DFACTOR                      AAXV0250
      WRITE(6,628) NGAM, NROINT, DRR                                    AAXV0251
      WRITE(6,629) RVDW, DVDW, ALVDW                                    AAXV0252
      WRITE(6,630) RDIA, DDIA, ALDIA                                    AAXV0253
      WRITE(6,622) IPRNT, IDWINTG, IDISC, IEDW                          AAXV0254
      IF (LWKBJ.NE.0) WRITE(6,623) STWKBJ                               AAXV0255
      WRITE(6,620) EPSIL                                                AAXV0256
      IF (IPRNT.GT.1) WRITE(6,649) KTOTAL,KSCR,KSR,KSI,KXK,KFK,KCRJG,   AAXV0257
     #                                            KETA,NCHM,NLSCR       AAXV0258
      WRITE(6,650) MAXC,NBYTE,KCR                                       AAXV0259
C                                                                       AAXV0260
      IF (EMIN.LT.0.0D0 .OR. ETOTAL.LT.0.0D0 .OR.                       AAXV0261
     +                     RMASS.LT.0.0D0 .OR. VIBQUAN.LT.0.0D0) THEN   AAXV0262
        WRITE(6,604)                                                    AAXV0263
        STOP 3                                                          AAXV0264
      END IF                                                            AAXV0265
C                                                                       AAXV0266
      IF (KCR.LE.0) THEN                                                AAXV0267
        WRITE(6,601) KCR, MAXC                                          AAXV0268
        STOP 2                                                          AAXV0269
      END IF                                                            AAXV0270
C                                                                       AAXV0271
      RETURN                                                            AAXV0272
C                                                                       AAXV0273
C INPUT-FORMATS                                                         AAXV0274
C                                                                       AAXV0275
  500 FORMAT(10A8)                                                      AAXV0276
  501 FORMAT(8D10.3)                                                    AAXV0277
C                                                                       AAXV0278
C OUTPUT-FORMATS                                                        AAXV0279
C                                                                       AAXV0280
  601 FORMAT(1H0,10(1H*),' PROBLEM IN CORE REQUEST : KCR = ',I3,        AAXV0281
     #                                      ' FOR ',I6,' WORDS')        AAXV0282
  604 FORMAT(1H0,11(1H*),' NO DEFAULT ALLOWED FOR EMIN - ETOTAL - ',    AAXV0283
     #                                      'RMASS - VIBQUAN ',11(1H*)) AAXV0284
  605 FORMAT( 1H1,2(/1X,129(1H*)),/,3H **,125X,2H**)                    AAXV0285
  606 FORMAT(3H **,40X,A5,4X,A6,4X,A5,4X,A6,4X,A7,40X,2H**)             AAXV0286
  607 FORMAT(3H **,20X,3(A6,4X),A7,12X,3(A6,4X),26X,2H**)               AAXV0287
  608 FORMAT(3H **,125X,2H**,2(/,1X,129(1H*)))                          AAXV0288
  611 FORMAT(1H0,20X,3H***,2X,10A8,2X,3H***)                            AAXV0289
  612 FORMAT(1H0/1H0,5X,'LIMITS ON VIBRATIONAL BASIS : NVMIN = ',I4,    AAXV0290
     #       ' NVMAX = ',I4)                                            AAXV0291
  613 FORMAT(1H0,5X,'ENERGY LIMITS : EMIN = ',D15.8,' (REDUCED :',D15.8 AAXV0292
     #       ,')','  ETOTAL = ',D15.8,' IN ',I4,' STEPS')               AAXV0293
  614 FORMAT(1H0,5X,'TOTAL ENERGY = ',D15.8,' (= ',D15.8,' REDUCED)')   AAXV0294
  615 FORMAT(1H0,5X,'TOTAL ANGULAR MOMENTUM LIMITS , JTMIN = ',I6,      AAXV0295
     #       ' JTMAX = ',I6,' JTSTEP = ',I5)                            AAXV0296
  617 FORMAT(1H0,5X,'ENERGY/DISTANCE SCALING : HBAR = ',D15.8,          AAXV0297
     #       ' RSCALE = ',D15.8,' VSCALE = ',D15.8)                     AAXV0298
  618 FORMAT(1H0,5X,'REDUCED MASS : ',D15.8,' (==> BRED = ',D15.8,')')  AAXV0299
  619 FORMAT(6X,'INTEGRATION RANGE : RBEG = ',D15.8,' (',D15.8,         AAXV0300
     #       ' REDUCED)',' REND = ',D15.8,' (',D15.8,' REDUCED)',       AAXV0301
     #      /,6X,20X,'DFACTOR = ',D10.3)                                AAXV0302
  620 FORMAT(1H0,5X,3(1H*),' THE PRECISION (EPSIL) IS ',D9.1,2X,3(1H*)) AAXV0303
  622 FORMAT(1H0,5X,'//OPTIONS  IPRNT :',I2,'  IDWINTG :',I2,'  IDISC :'AAXV0304
     #       I3,'   IEDW :',I3)                                         AAXV0305
  623 FORMAT(1H0,5X,'// A WKBJ-PROPAGATOR IS USED AS SOON AS ',         AAXV0306
     #       ' P',1H','/2*P**2 < ',D12.4,3X,2H//)                       AAXV0307
  624 FORMAT(6X,'SCALE-FACTOR FOR TOTAL AND PARTIAL J CROSSECTIONS'     AAXV0308
     #       ,' RESPECTIVELY : ',D15.8,' / ',D15.8,/,6X,                AAXV0309
     #       '( A ZERO YIELDS THE DIMENSIONLESS CROSSECTION)')          AAXV0310
  628 FORMAT(1H0,5X,'CALCULATIONS OVER  ',I6,'  SCATTERING ',           AAXV0311
     &           ' ANGLES ','    NUMBER OF GRID POINTS = ',I4,          AAXV0312
     &                           ' IN STEPS OF ',D14.6)                 AAXV0313
  629 FORMAT(1H0,5X,'VDW POTENTIAL: R0= ',D15.8,'  D= ',D15.8,          AAXV0314
     &                             '  ALPHA =',D15.8)                   AAXV0315
  630 FORMAT(1H0,5X,'DIA POTENTIAL: R0= ',D15.8,'  D= ',D15.8,          AAXV0316
     &                                        '  ALPHA= ',D15.8)        AAXV0317
  636 FORMAT(3H **,43X,A5,5X,A6,5X,A6,5X,A7,43X,2H**)                   AAXV0318
  637 FORMAT(3H **,25X,2(A6,4X),A7,8X,3(A6,4X),33X,2H**)                AAXV0319
  649 FORMAT(1H0,5X,8HKTOTAL :,I6,3X,8HKSCR   :,I6,3X,8HKSR    :,I6,3X, AAXV0320
     #       8HKSI    :,I6,/,6X,8HKXK    :,I6,3X,8HKFK    :,I6,         AAXV0321
     #       3X,8HKCRJG  :,I6,3X,8HKETA   :,I6,/,6X,                    AAXV0322
     #       8HNCHM   :,I6,3X,8HNLSCR  :,I6)                            AAXV0323
  650 FORMAT(1H0/1H0,5X,10(1H-),' THIS RUN REQUIRES ',I6,               AAXV0324
     #       ' WORDS OF BLANK COMMON (',I2,' BYTES / KCR = ',I6,' )',   AAXV0325
     #       10(1H-),////,1X,129(1H*))                                  AAXV0326
C                                                                       AAXV0327
      END                                                               AAXV0328
C                                                                       AAXV0329
      SUBROUTINE EDW(SCR,SREAL,SIMAG,ETA,NJ,NROW,FLUX,XK,FK,UCF,NBB,    AAXV0330
     #                                CRJG,CRJ,TOTAL,RR,POTRX,DPOTRX)   AAXV0331
C********************************************************************   AAXV0332
C                                                                       AAXV0333
C DRIVER PROGRAM FOR CALCULATION OF CROSS SECTIONS                      AAXV0334
C THREE FORTRAN CHANNELS ARE USED:                                      AAXV0335
C       10 - POTENTIAL MATRIX                                           AAXV0336
C       11 - ITS DERIVATIVES                                            AAXV0337
C       12 - PROPAGATION GRID POINTS                                    AAXV0338
C                                                                       AAXV0339
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0340
      INTEGER*4 POINTR,CUE,PT                                           AAXV0341
      EXTERNAL D01BAZ                                                   AAXV0342
C                                                                       AAXV0343
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV0344
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV0345
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV0346
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV0347
     #                                        NVMIN,NVMAX               AAXV0348
      COMMON/SURF / IPNT,NPOINT,NVMXP1                                  AAXV0349
      COMMON/ANG  / DRR,NGAM,NROINT,IOLD                                AAXV0350
C                                                                       AAXV0351
      PARAMETER(NDIM=96,N2DIM=4)                                        AAXV0352
C                                                                       AAXV0353
      DIMENSION SCR(1),SREAL(1),SIMAG(1),ETA(1),NJ(1),NROW(1),          AAXV0354
     #          FLUX(1),XK(1),FK(1),UCF(1),NBB(1),CRJG(1),CRJ(1),       AAXV0355
     *          TOTAL(1),XVALUE(NDIM),VIBM(NDIM,N2DIM),POTRX(1),        AAXV0356
     *          J(1),RR(1),DPOTRX(1),WEIGHT(20),ABSCI(20)               AAXV0357
C                                                                       AAXV0358
      DATA ZERO/0.0D0/                                                  AAXV0359
      DATA POINTR,PT/'STRT','END.'/                                     AAXV0360
C                                                                       AAXV0361
      IF (N2DIM.NE.NVMXP1) THEN                                         AAXV0362
        WRITE(6,600)                                                    AAXV0363
        RETURN                                                          AAXV0364
      END IF                                                            AAXV0365
C                                                                       AAXV0366
      IPNT   = 1                                                        AAXV0367
      ALIM   =-1.0D0                                                    AAXV0368
      BLIM   = 1.0D0                                                    AAXV0369
      ITYPE  = 0                                                        AAXV0370
      IFAIL  = 0                                                        AAXV0371
      NPOINT = NROINT                                                   AAXV0372
      NMAX   = NPOINT*NVMXP1*(NVMXP1+1)/2                               AAXV0373
C  FIND WEIGHTS AND ABSCISSAE FOR GAUSSIAN QUADRATURE                   AAXV0374
      CALL D01BBF(D01BAZ,ALIM,BLIM,ITYPE,NGAM,WEIGHT,ABSCI,IFAIL)       AAXV0375
C                                                                       AAXV0376
      IF (IOLD.NE.1) THEN                                               AAXV0377
C  EVALUATE POTENTIAL MATRIX AND DERIVATIVES                            AAXV0378
C  FIRST FIND VIBRATIONAL WAVE FUNCTIONS                                AAXV0379
        CALL VIBR(XVALUE,VIBM)                                          AAXV0380
        OPEN (UNIT=10,FILE='POT.DAT',STATUS='NEW',FORM='UNFORMATTED')   AAXV0381
        OPEN (UNIT=11,FILE='DPOT.DAT',STATUS='NEW',FORM='UNFORMATTED')  AAXV0382
        OPEN (UNIT=12,FILE='GRID.DAT',STATUS='NEW',FORM='UNFORMATTED')  AAXV0383
        REWIND 10                                                       AAXV0384
        REWIND 11                                                       AAXV0385
        REWIND 12                                                       AAXV0386
C                                                                       AAXV0387
        RR(1)  = RMIN                                                   AAXV0388
        DO 10 I=1,NPOINT                                                AAXV0389
        RR(I+1)= RR(I)+DRR                                              AAXV0390
   10   CONTINUE                                                        AAXV0391
        WRITE(12) (RR(I), I=1,NPOINT)                                   AAXV0392
C                                                                       AAXV0393
        DO 25 I=1,NGAM,2                                                AAXV0394
        GAMMA  = DABS(DACOS(ABSCI(I)))                                  AAXV0395
        WRITE(10) POINTR                                                AAXV0396
        WRITE(10) GAMMA                                                 AAXV0397
        WRITE(11) POINTR                                                AAXV0398
        WRITE(11) GAMMA                                                 AAXV0399
C FIND POTENTIAL MATRIX AND DERIVATIVES                                 AAXV0400
        CALL POTMAT(RR,GAMMA,XVALUE,POTRX,DPOTRX,NPOINT,VIBM)           AAXV0401
        WRITE(10) (POTRX(IJK), IJK=1,NMAX)                              AAXV0402
        WRITE(11) (DPOTRX(IJK), IJK=1,NMAX)                             AAXV0403
   25   CONTINUE                                                        AAXV0404
        WRITE(10) PT                                                    AAXV0405
C WRITING OF POTENTIAL MATRIX TO DISK IS NOW COMPLETE                   AAXV0406
        CLOSE (UNIT=10)                                                 AAXV0407
        CLOSE (UNIT=11)                                                 AAXV0408
        CLOSE (UNIT=12)                                                 AAXV0409
C                                                                       AAXV0410
      END IF                                                            AAXV0411
C                                                                       AAXV0412
C IF POTENTIAL MATRIX EXISTS, CONTINUE                                  AAXV0413
C                                                                       AAXV0414
        OPEN (UNIT=10,FILE='POT.DAT',STATUS='OLD',FORM='UNFORMATTED')   AAXV0415
        OPEN (UNIT=11,FILE='DPOT.DAT',STATUS='OLD',FORM='UNFORMATTED')  AAXV0416
        OPEN (UNIT=12,FILE='GRID.DAT',STATUS='OLD',FORM='UNFORMATTED')  AAXV0417
C                                                                       AAXV0418
      REWIND 12                                                         AAXV0419
      IF (IDISC.NE.0)              REWIND IDISC                         AAXV0420
      READ(12) (RR(I), I=1,NPOINT)                                      AAXV0421
      NTSTEP = (JJTMAX-JJTMIN)/JTSTEP + 1                               AAXV0422
C                                                                       AAXV0423
C LOOP OVER TOTAL ENERGY                                                AAXV0424
C                                                                       AAXV0425
      DO 450 NEINT=1,NETOT                                              AAXV0426
      ERED = ETOT + DFLOAT(NEINT-1)*ESTEP                               AAXV0427
C                                                                       AAXV0428
C CALCULATE CHANNEL INTERNAL ENERGIES, WAVE NUMBERS, FLUX FACTORS       AAXV0429
C                                                                       AAXV0430
      CALL ELEVEL(BRED,ERED,VIBQUAN,NVMIN,NVMAX,                        AAXV0431
     &                                 NOPEN,KJ,SCR,XK,FK,FLUX)         AAXV0432
C                                                                       AAXV0433
      KJMAX = KJ*(KJ+1)/2                                               AAXV0434
C                                                                       AAXV0435
C INITIALISE TOTAL CROSS SECTIONS                                       AAXV0436
C                                                                       AAXV0437
      DO 100 K = 1,KJMAX                                                AAXV0438
  100 TOTAL(K) = ZERO                                                   AAXV0439
C INITIALISE TOTAL ANGULAR MOMENTUM INTERVAL INDEX                      AAXV0440
      JN   = 0                                                          AAXV0441
      XCRJ = 0.0D0                                                      AAXV0442
C LOOP OVER GOOD QUANTUM NUMBERS                                        AAXV0443
      DO 400 JJTOT = JJTMIN,JJTMAX,JTSTEP                               AAXV0444
      REWIND 10                                                         AAXV0445
      REWIND 11                                                         AAXV0446
      JTOT = JJTOT-1                                                    AAXV0447
C INITIALISE J CROSS SECTIONS                                           AAXV0448
      DO 110 K = 1,KJMAX                                                AAXV0449
  110 CRJ(K)   = ZERO                                                   AAXV0450
      CONTINUE                                                          AAXV0451
C                                                                       AAXV0452
C LOOP OVER SCATTERING ANGLES                                           AAXV0453
      DO 300 KK = 1,NGAM,2                                              AAXV0454
C READ POTENTIAL AND ITS DERIVATIVES FROM FILES                         AAXV0455
      IPNT  = 1                                                         AAXV0456
      GAMMA = DABS(DACOS(ABSCI(KK)))                                    AAXV0457
  152 READ(10) CUE                                                      AAXV0458
      IF(CUE.EQ.PT) THEN                                                AAXV0459
        WRITE (6,640)                                                   AAXV0460
        STOP 55                                                         AAXV0461
      END IF                                                            AAXV0462
      IF(CUE.NE.POINTR)                       GOTO 152                  AAXV0463
      READ (10) ANGLE                                                   AAXV0464
      IF (DABS(ANGLE-GAMMA).GT.0.0001D0)      GOTO 152                  AAXV0465
      READ(10) (POTRX(I), I = 1,NMAX)                                   AAXV0466
C                                                                       AAXV0467
  155 READ(11) CUE                                                      AAXV0468
      IF(CUE.EQ.PT)                     STOP 55                         AAXV0469
      IF(CUE.NE.POINTR)                       GOTO 155                  AAXV0470
      READ (11) ANGLE                                                   AAXV0471
      IF (DABS(ANGLE-GAMMA) .GT. 0.0001D0)    GOTO 155                  AAXV0472
      READ(11) (DPOTRX(I), I = 1,NMAX)                                  AAXV0473
C                                                                       AAXV0474
      N1  = 0                                                           AAXV0475
      NNJ = 0                                                           AAXV0476
      DO 180 IV = NVMIN,NVMAX                                           AAXV0477
      NNJ = NNJ+1                                                       AAXV0478
      IF(NNJ.GT.KJ)                     GOTO 185                        AAXV0479
      N1      = N1+1                                                    AAXV0480
      J(N1)   = IV + 1                                                  AAXV0481
      UCF(N1) = DFLOAT(JTOT)                                            AAXV0482
      NBB(N1) = NNJ                                                     AAXV0483
      NROW(N1) = NNJ                                                    AAXV0484
  180 CONTINUE                                                          AAXV0485
  185 NCHAN   = N1                                                      AAXV0486
C                                                                       AAXV0487
      IF (NCHAN.NE.0) THEN                                              AAXV0488
        DO 190 N1 = 1,NCHAN                                             AAXV0489
        UCF(N1) = UCF(N1)*(UCF(N1)+1.D0)                                AAXV0490
  190   CONTINUE                                                        AAXV0491
        IF (IPRNT.GE.2) THEN                                            AAXV0492
          WRITE(6,656) JTOT, J(1), J(NCHAN), NCHAN                      AAXV0493
          WRITE(6,660)                                                  AAXV0494
          DO 195 N1 = 1,NCHAN                                           AAXV0495
          XL = (DSQRT(1.D0 +4.D0*UCF(N1))-1.D0)/2.D0                    AAXV0496
          WRITE(6,662) N1,J(N1),XL,UCF(N1)                              AAXV0497
  195     CONTINUE                                                      AAXV0498
        END IF                                                          AAXV0499
      END IF                                                            AAXV0500
C                                                                       AAXV0501
      IF (NCHAN.NE.0) THEN                                              AAXV0502
C                                                                       AAXV0503
C CALCULATE T-MATRIX                                                    AAXV0504
C                                                                       AAXV0505
        CALL TMAT(SREAL,SIMAG,ETA,NJ,NROW,XK,FK,SCR,UCF,NBB,RR,         AAXV0506
     &                                                  POTRX,DPOTRX)   AAXV0507
C                                                                       AAXV0508
C CALCULATE CROSS-SECTIONS                                              AAXV0509
C                                                                       AAXV0510
C INITIALIZE GAMMA CROSS-SECTIONS                                       AAXV0511
C                                                                       AAXV0512
        DO 210 I = 1,KJMAX                                              AAXV0513
  210   CRJG(I)  = ZERO                                                 AAXV0514
        CALL CROSS(NCHAN,NROW,SREAL,SIMAG,CRJG)                         AAXV0515
        DO 220 I=1,KJMAX                                                AAXV0516
  220   CRJ(I)   = CRJ(I) + CRJG(I)*WEIGHT(KK)*2.0D0                    AAXV0517
C                                                                       AAXV0518
        IF (IPRNT.GT.1) THEN                                            AAXV0519
          WRITE(6,609) JTOT,GAMMA                                       AAXV0520
          CALL WRITEC(1,KJ,1,KJ,CONVP,FLUX,XK,CRJG)                     AAXV0521
        END IF                                                          AAXV0522
      END IF                                                            AAXV0523
  300 CONTINUE                                                          AAXV0524
C                                                                       AAXV0525
C  DETERMINE WEIGHT FOR THIS POINT, USING VARIOUS QUADRATURE SCHEMES    AAXV0526
C                                                                       AAXV0527
      DJTOT    = QUAD(JN,NTSTEP,JTSTEP)                                 AAXV0528
      DO 390 I = 1,KJMAX                                                AAXV0529
  390 TOTAL(I) = TOTAL(I) + CRJ(I)*DJTOT                                AAXV0530
C                                                                       AAXV0531
C CLOSE LOOP ON TOTAL J                                                 AAXV0532
C                                                                       AAXV0533
      IF (IPRNT.GT.0) THEN                                              AAXV0534
        WRITE(6,610) JTOT,ERED                                          AAXV0535
        CALL WRITEC(1,KJ,1,KJ,CONVP,FLUX,XK,CRJ)                        AAXV0536
      END IF                                                            AAXV0537
      IF(CRJ(2).GT.XCRJ)           XCRJ = CRJ(2)                        AAXV0538
      IF(CRJ(2).LT.XCRJ*0.05D0)            GOTO 425                     AAXV0539
C                                                                       AAXV0540
  400 CONTINUE                                                          AAXV0541
  425 CONTINUE                                                          AAXV0542
      XPRO = 100.0D0*CRJ(2)/TOTAL(2)                                    AAXV0543
      CCR  = CRJ(2)*FLUX(1)                                             AAXV0544
      WRITE(6,620) JTOT, CCR, XPRO                                      AAXV0545
C                                                                       AAXV0546
      JTMIN = JJTMIN - 1                                                AAXV0547
      JTMAX = JJTMAX - 1                                                AAXV0548
      WRITE(6,611) BRED,ERED,JTMIN,JTOT,JTSTEP                          AAXV0549
      CALL WRITEC(1,KJ,1,KJ,CONV,FLUX,XK,TOTAL)                         AAXV0550
      WRITE(6,612)                                                      AAXV0551
C                                                                       AAXV0552
C CLOSE LOOP ON TOTAL ENERGY                                            AAXV0553
C                                                                       AAXV0554
  450 CONTINUE                                                          AAXV0555
C                                                                       AAXV0556
  600   FORMAT(1H0,10(1H*),' PARAMETER N2DIM IS NOT EQUAL TO ',         AAXV0557
     ,                                          'NVMAX+1 ',10(1H*))     AAXV0558
  609   FORMAT(///,15X,'---- INDIVIDUAL GAMMA CROSSECTIONS  ----',      AAXV0559
     #                      20X,'JTOT ',I6,10X,'GAMMA=  ',F10.5)        AAXV0560
  610   FORMAT(//,10X,'---- INDIVIDUAL J CROSSECTIONS ----',16X,        AAXV0561
     #       '::::  JTOT ',I4,' ::::::  ERED  ',D12.6,' ::::')          AAXV0562
  611   FORMAT(//,1X,130(1H-),//,15X,'====== TOTAL CROSSECTIONS ======',AAXV0563
     #       //,10X,'BRED ',D12.5,5X,'ERED ',D12.5,5X,                  AAXV0564
     #       5X,'JTOT FROM ',I6,'   TO ',I6,'  WITH STEP ',             AAXV0565
     #       I6,//,1X)                                                  AAXV0566
  612   FORMAT(1H0,/,1X,130(1H-),//)                                    AAXV0567
  620   FORMAT(//,10X,'FINAL CONTRIBUTION FROM JTOT ',I6,' WAS ',       AAXV0568
     &             D14.7,' REPRESENTING ',D10.4, '% OF TOTAL',/)        AAXV0569
  640   FORMAT(1H ,5X,'MISMATCH BETWEEN PARAMETER NGAM AND ',           AAXV0570
     &         'ITS VALUE ON DISK')                                     AAXV0571
  656   FORMAT(1H0,//,'     JTOT= ',I5, ' MIN LEVEL= ',I3,              AAXV0572
     *          '     MAX LEVEL= ',I3,4X,' NO. OF CHANNELS= ',I3)       AAXV0573
  660   FORMAT('O     N1    LEVEL         L  ',8X,'UCF',/)              AAXV0574
  662   FORMAT(1H ,2X,2(I4,4X),3X,F6.0,F12.1)                           AAXV0575
C                                                                       AAXV0576
      RETURN                                                            AAXV0577
      END                                                               AAXV0578
C                                                                       AAXV0579
      SUBROUTINE ATOS(SREAL,SIMAG,ETA,NCHAN)                            AAXV0580
C*******************************************************************    AAXV0581
C                                                                       AAXV0582
C COMPUTE FROM THE ACTION MATRIX (0,SIMAG)                              AAXV0583
C             THE DISTORTED-WAVE S-MATRIX (SREAL,SIMAG)                 AAXV0584
C                                                                       AAXV0585
C NCHAN  : NUMBER OF CHANNELS                                           AAXV0586
C ETA    : PHASE SHIFTS                                                 AAXV0587
C                                                                       AAXV0588
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0589
      DIMENSION SREAL(1),SIMAG(1),ETA(1)                                AAXV0590
C                                                                       AAXV0591
      N = 0                                                             AAXV0592
      DO 100 N1 = 1,NCHAN                                               AAXV0593
      DO 50  N2 = 1,N1                                                  AAXV0594
      N = N + 1                                                         AAXV0595
      ETAS     = ETA(N1) + ETA(N2)                                      AAXV0596
      SREAL(N) =-SIMAG(N)*DSIN(ETAS)                                    AAXV0597
      SIMAG(N) = SIMAG(N)*DCOS(ETAS)                                    AAXV0598
   50 CONTINUE                                                          AAXV0599
C                                                                       AAXV0600
      SREAL(N) = SREAL(N) + DCOS(ETAS)                                  AAXV0601
      SIMAG(N) = SIMAG(N) + DSIN(ETAS)                                  AAXV0602
  100 CONTINUE                                                          AAXV0603
C                                                                       AAXV0604
      RETURN                                                            AAXV0605
      END                                                               AAXV0606
C                                                                       AAXV0607
      SUBROUTINE CAN(X,C0,C1,C2,C3)                                     AAXV0608
C*******************************************************************    AAXV0609
C                                                                       AAXV0610
C CN(X) IS (1/2)*INTEGRAL OVER T FROM -1 TO +1 OF                       AAXV0611
C                            (T**N)*SIN(T*X)   IF N IS ODD              AAXV0612
C                            (T**N)*COS(T*X)   IF N IS EVEN             AAXV0613
C                                                                       AAXV0614
C CALLED BY QQINT                                                       AAXV0615
C                                                                       AAXV0616
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0617
      DATA XMIN /1.0D-5/                                                AAXV0618
      DATA D3  /0.333333333333D0/                                       AAXV0619
      DATA D6  /0.166666666667D0/                                       AAXV0620
      DATA D30 /3.333333333333D-2/                                      AAXV0621
C                                                                       AAXV0622
      IF(DABS(X) .LE. XMIN) THEN                                        AAXV0623
        X2 = X*X                                                        AAXV0624
        X3 = X*X2                                                       AAXV0625
        C0 = 1.0D0-X2*D6                                                AAXV0626
        C1 = X*D3-X3*D30                                                AAXV0627
        C2 = 1.0D0/3.0D0-0.1D0*X2                                       AAXV0628
        C3 = X*(0.25D0-X2/42.0D0)                                       AAXV0629
      ELSE                                                              AAXV0630
        S  = DSIN(X)                                                    AAXV0631
        C  = DCOS(X)                                                    AAXV0632
        XINV = 1.0D0/X                                                  AAXV0633
        C0 = S*XINV                                                     AAXV0634
        C1 = (C0-C)*XINV                                                AAXV0635
        C2 = (S-C1-C1)*XINV                                             AAXV0636
        C3 = (3.0D0*C2-C)*XINV                                          AAXV0637
      END IF                                                            AAXV0638
C                                                                       AAXV0639
      RETURN                                                            AAXV0640
      END                                                               AAXV0641
C                                                                       AAXV0642
      SUBROUTINE CROSS(MAX,MROW,TREAL,TIMAG,CRJ)                        AAXV0643
C*******************************************************************    AAXV0644
C                                                                       AAXV0645
C COMPUTE CROSSECTIONS FROM T-MATRIX                                    AAXV0646
C ONLY THE PART , WHICH IS THE SAME FOR NJ=>NJ' AND NJ'=>NJ             AAXV0647
C IS COMPUTED HERE; THE REST IS DONE IN WRITEC                          AAXV0648
C                                                                       AAXV0649
          IMPLICIT REAL*8(A-H,O-Z)                                      AAXV0650
          DIMENSION MROW(1),TREAL(1),TIMAG(1),CRJ(1)                    AAXV0651
C                                                                       AAXV0652
          COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ        AAXV0653
C                                                                       AAXV0654
          FL1 = 2*JTOT+1                                                AAXV0655
          M  = 0                                                        AAXV0656
          DO 80 NJ1=1,MAX                                               AAXV0657
          J1 = MROW(NJ1)                                                AAXV0658
          DO 80 NJ2=1,NJ1                                               AAXV0659
          J2 = MROW(NJ2)                                                AAXV0660
          M  = M + 1                                                    AAXV0661
          IF(J2.GT.J1) THEN                                             AAXV0662
            NNJ = (J2*J2-J2)/2 + J1                                     AAXV0663
          ELSE                                                          AAXV0664
            NNJ = (J1*J1-J1)/2 + J2                                     AAXV0665
          END IF                                                        AAXV0666
          CR = (TREAL(M)**2+TIMAG(M)**2)*FL1/2.0D0                      AAXV0667
          CRJ(NNJ) = CRJ(NNJ)+CR                                        AAXV0668
   80     CONTINUE                                                      AAXV0669
C                                                                       AAXV0670
      RETURN                                                            AAXV0671
      END                                                               AAXV0672
C                                                                       AAXV0673
      SUBROUTINE ELEVEL(BRED,ERED,VIBQUAN,NVMIN,NVMAX,                  AAXV0674
     &                                  NOPEN,KJ,XKK,XK,FK,FLUX)        AAXV0675
C*******************************************************************    AAXV0676
C                                                                       AAXV0677
C CALCULATES CHANNEL INTERNAL ENERGIES, WAVENUMBERS, FLUX               AAXV0678
C NORMALIZATION FACTORS AND NUMBER OF OPEN CHANNELS                     AAXV0679
C                                                                       AAXV0680
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0681
      DIMENSION XK(1),FK(1),XKK(1),FLUX(1)                              AAXV0682
C                                                                       AAXV0683
      DATA ZERO, ONE /0.0D0,1.0D0/                                      AAXV0684
      DATA PI /3.14159265358979D0/                                      AAXV0685
C                                                                       AAXV0686
      WRITE(6,612)                                                      AAXV0687
      NNVMIN = NVMIN+1                                                  AAXV0688
      NNVMAX = NVMAX+1                                                  AAXV0689
      NNJ   = 0                                                         AAXV0690
      NOPEN = 0                                                         AAXV0691
      DO 120 NNV = NNVMIN,NNVMAX                                        AAXV0692
      NV1   = NNV-1                                                     AAXV0693
      NNJ   = NNJ+1                                                     AAXV0694
      J1    = JJ-1                                                      AAXV0695
      EINT  = VIBQUAN*(NV1+0.5D0)-0.6127*4.554D-6*NV1**2                AAXV0696
      XKK1  = BRED*(ERED-EINT)                                          AAXV0697
      XKK(NNJ) = XKK1                                                   AAXV0698
      XK1   = DSQRT(DABS(XKK1))                                         AAXV0699
      XK(NNJ) = XK1                                                     AAXV0700
      IF(XKK1 .GT. ZERO)    NOPEN = NOPEN+1                             AAXV0701
      FK(NNJ)   = ONE/DSQRT(XK1)                                        AAXV0702
      FLUX(NNJ) = PI/XKK1                                               AAXV0703
      XKINET    = XKK1/BRED                                             AAXV0704
      WRITE(6,605) NNJ,NV1,EINT,XKK1,XKINET                             AAXV0705
  120 CONTINUE                                                          AAXV0706
      KJ = NNJ                                                          AAXV0707
C CHECK THAT THERE ARE OPEN CHANNELS                                    AAXV0708
      WRITE(6,604) NOPEN,NNJ                                            AAXV0709
      IF(NOPEN .LT. KJ) KJ = NOPEN                                      AAXV0710
C                                                                       AAXV0711
  604 FORMAT(1H0,/,5X,'NUMBER OF OPEN CHANNELS=',I3,' TOTAL NUMBER=',   AAXV0712
     &        I3/)                                                      AAXV0713
  605 FORMAT(1H ,2I6,3(2X,D15.7))                                       AAXV0714
  612 FORMAT(1H0,///,25X,'--ELEVEL--',//,5X,'NJ',4X,'NV',               AAXV0715
     #           6X,'VIBRATION',8X,'W. VECTOR',7X,'  KINETIC')          AAXV0716
      RETURN                                                            AAXV0717
      END                                                               AAXV0718
C                                                                       AAXV0719
      REAL FUNCTION QRB*8(L,X)                                          AAXV0720
C*******************************************************************    AAXV0721
C                                                                       AAXV0722
C UNIFORM APPROX TO PHASE OF THE REGULAR RICCATI-BESSEL FUNCTION        AAXV0723
C                                                                       AAXV0724
      IMPLICIT REAL*8(A-H,O-Z)                                          AAXV0725
      REAL*8 L                                                          AAXV0726
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV0727
      PI4 = DATAN(1.0D0)                                                AAXV0728
      IF(DABS(L).LT.1.0D-20) THEN                                       AAXV0729
        QRB = X                                                         AAXV0730
      ELSE                                                              AAXV0731
        XL = L+0.5D0                                                    AAXV0732
        IF(X .LT. XL) THEN                                              AAXV0733
C FORBIDDEN REGION                                                      AAXV0734
          XM = XL + 0.5D0                                               AAXV0735
          QRB= (X+(XM*PI4-XL)*(X/XL)**2)/XM                             AAXV0736
        ELSE                                                            AAXV0737
          ROOT = DSQRT((X/XL)**2-1.0D0)                                 AAXV0738
          QRB  = PI4+XL*(ROOT-DATAN(ROOT))/DSQRT(BRED)                  AAXV0739
        END IF                                                          AAXV0740
      END IF                                                            AAXV0741
      RETURN                                                            AAXV0742
      END                                                               AAXV0743
C                                                                       AAXV0744
      SUBROUTINE RENORM(NCHAN,APOS,ANEG,AMAT)                           AAXV0745
C*******************************************************************    AAXV0746
C                                                                       AAXV0747
C RENORMALIZE AMPLITUDES OF WAVEFUNCTIONS AND THE A-MATRIX              AAXV0748
C                                                                       AAXV0749
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0750
      DIMENSION APOS(1),ANEG(1),AMAT(1)                                 AAXV0751
      COMMON /OPT  / IPRNT,IDWINTG,IDISC                                AAXV0752
C                                                                       AAXV0753
      DATA APZ,AAPZ/1.0D-20,1.0D-30/                                    AAXV0754
C                                                                       AAXV0755
      IF(IPRNT.GT.1)             WRITE(6,607)                           AAXV0756
      K = 0                                                             AAXV0757
      DO 227 K1 = 1,NCHAN                                               AAXV0758
      APOS(K1)  = APZ*APOS(K1)                                          AAXV0759
      ANEG(K1)  = APZ*ANEG(K1)                                          AAXV0760
      IF (IPRNT.GT.2)            WRITE(6,600) K1,APOS(K1),ANEG(K1)      AAXV0761
      DO 227 K2 = 1,K1                                                  AAXV0762
      K = K+1                                                           AAXV0763
  227 AMAT(K) = AAPZ*AMAT(K)                                            AAXV0764
C                                                                       AAXV0765
      IF (IPRNT.LT.3)            RETURN                                 AAXV0766
      WRITE(6,601)                                                      AAXV0767
      CALL WRITER(AMAT,NCHAN)                                           AAXV0768
C                                                                       AAXV0769
  600 FORMAT(5X,' N/AMPLITUDES (POS/NEG) : ',I5,2D20.5)                 AAXV0770
  601 FORMAT(1H0,5X,'  THE RENORMALIZED A-MATRIX ')                     AAXV0771
  607 FORMAT(1H0,5(1H*),' THE WAVEFUNCTIONS AND THE A-MATRIX ',         AAXV0772
     #       'ARE RENORMALIZED  (*E-20/E-40 RESPECTIVELY)')             AAXV0773
C                                                                       AAXV0774
      RETURN                                                            AAXV0775
      END                                                               AAXV0776
C                                                                       AAXV0777
      SUBROUTINE TMAT(SREAL,SIMAG,ETA,NJ,NROW,XK,FK,SCR,UCF,NBB,RR,     AAXV0778
     &                                             POTRX,DPOTRX)        AAXV0779
C*******************************************************************    AAXV0780
C                                                                       AAXV0781
C CALCULATE T-MATRIX                                                    AAXV0782
C                                                                       AAXV0783
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0784
C                                                                       AAXV0785
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV0786
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV0787
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV0788
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV0789
     #              NVMIN,NVMAX                                         AAXV0790
C                                                                       AAXV0791
      DATA ZERO/0.0D0/                                                  AAXV0792
C                                                                       AAXV0793
      DIMENSION SREAL(1),SIMAG(1),ETA(1),NJ(1),NROW(1),XK(1),SCR(1),    AAXV0794
     #          FK(1),UCF(1),NBB(1),RR(1),POTRX(1),DPOTRX(1)            AAXV0795
C                                                                       AAXV0796
C     INITIALIZE S-MATRICES AND PHASE-SHIFTS                            AAXV0797
C                                                                       AAXV0798
      NLIN      = NCHAN*(NCHAN+1)/2                                     AAXV0799
      DO 200 N  = 1,NLIN                                                AAXV0800
  200 SIMAG(N)  = ZERO                                                  AAXV0801
      DO 205 N1 = 1,NCHAN                                               AAXV0802
  205 ETA(N1)   = ZERO                                                  AAXV0803
C                                                                       AAXV0804
C     K = LENGTH+1 OF SCRATCH-ARRAYS IN ACTION  I.E.                    AAXV0805
C     RSTART, SP, SN, CP, CN, ANEG, QNEG, APOS, QPOS, WN, WP,           AAXV0806
C     BB, AND TP ARE ALL LOCATED IN SCR                                 AAXV0807
C     K IS ALSO USED FOR THE CORE-PARTITIONING IN EXPA                  AAXV0808
C                                                                       AAXV0809
      K = NCHAN + 1                                                     AAXV0810
C                                                                       AAXV0811
      CALL ACTION(ETA,SIMAG,XK,FK,NROW,UCF,NBB,                         AAXV0812
     #            SCR,SCR(K),SCR(2*K),SCR(3*K),SCR(4*K),SCR(5*K),       AAXV0813
     #            SCR(6*K),SCR(7*K),SCR(8*K),SCR(9*K),SCR(10*K),        AAXV0814
     #            SCR(11*K),SCR(12*K),NVIB,MVIB,RR,POTRX,DPOTRX)        AAXV0815
C                                                                       AAXV0816
      IF(IPRNT .GE. 2) THEN                                             AAXV0817
C PRINT A-MATRIX                                                        AAXV0818
        WRITE(6,601) JTOT, ERED                                         AAXV0819
        CALL WRITER(SIMAG,NCHAN)                                        AAXV0820
      END IF                                                            AAXV0821
C                                                                       AAXV0822
C CALCULATE S-MATRIX                                                    AAXV0823
C                                                                       AAXV0824
C      EXPONENTIAL DISTORTED WAVE                                       AAXV0825
C                                                                       AAXV0826
       IF (IABS(IEDW).NE.2) THEN                                        AAXV0827
         CALL EXPA(SREAL,SIMAG,ETA,NCHAN,SCR,                           AAXV0828
     &            SCR(NCHAN*K),SCR((NCHAN+5)*K),SCR(NCHAN*K),           AAXV0829
     &                      SCR((NCHAN+1)*K),SCR((NCHAN+10)*K))         AAXV0830
C                                                                       AAXV0831
C      DISTORTED WAVE                                                   AAXV0832
C                                                                       AAXV0833
       ELSE                                                             AAXV0834
         CALL ATOS(SREAL,SIMAG,ETA,NCHAN)                               AAXV0835
       END IF                                                           AAXV0836
C                                                                       AAXV0837
       IF(IPRNT.GE.2) THEN                                              AAXV0838
C                                                                       AAXV0839
C  PRINT PHASESHIFTS AND S-MATRICES                                     AAXV0840
C                                                                       AAXV0841
         WRITE(6,603) (ETA(K),K=1,NCHAN)                                AAXV0842
         WRITE(6,604)                                                   AAXV0843
         CALL WRITER(SREAL,NCHAN)                                       AAXV0844
         WRITE(6,605)                                                   AAXV0845
         CALL WRITER(SIMAG,NCHAN)                                       AAXV0846
       END IF                                                           AAXV0847
C                                                                       AAXV0848
       IF(IDISC.NE.0) THEN                                              AAXV0849
C WRITE S-MATRICES TO FILE IDISC                                        AAXV0850
         NLIN = NCHAN*(NCHAN+1)/2                                       AAXV0851
         WRITE(6,606) NLIN, JTOT                                        AAXV0852
         WRITE(IDISC) JTOT, NLIN                                        AAXV0853
         WRITE(IDISC) (SREAL(K),K=1,NLIN)                               AAXV0854
         WRITE(IDISC) (SIMAG(K),K=1,NLIN)                               AAXV0855
       END IF                                                           AAXV0856
C                                                                       AAXV0857
C CALCULATE T-MATRIX FROM S-MATRIX                                      AAXV0858
C                                                                       AAXV0859
      DO 300 N1 = 1,NCHAN                                               AAXV0860
      NDIAG     = N1*(N1+1)/2                                           AAXV0861
  300 SREAL(NDIAG) = SREAL(NDIAG)-1.0D0                                 AAXV0862
C                                                                       AAXV0863
  601 FORMAT(1H0,5X,7(1H-),' A - MATRIX ',7(1H-),10X,                   AAXV0864
     #       8X,' JTOT ',I6,' ERED ',D15.7)                             AAXV0865
  603 FORMAT(1H0,/,6X,7(1H-),' // PHASE-SHIFTS // ',7(1H-),//,          AAXV0866
     #       (8X,8F15.7))                                               AAXV0867
  604 FORMAT(1H0,/,6X,7(1H-),'    REAL   PART OF S-MATRIX ',7(1H-))     AAXV0868
  605 FORMAT(1H0,/,6X,7(1H-),' IMAGINARY PART OF S-MATRIX ',7(1H-))     AAXV0869
  606 FORMAT(1H0,'MATRICES WRITTEN TO FILE / NELEMNTS,JTOT : ',2I10)    AAXV0870
C                                                                       AAXV0871
      RETURN                                                            AAXV0872
      END                                                               AAXV0873
C                                                                       AAXV0874
      SUBROUTINE EXPA(SREAL,SIMAG,ETA,NCHAN,QMAT,SCR,ADIAG,CA,SA,VECTR) AAXV0875
C********************************************************************** AAXV0876
C                                                                       AAXV0877
C EXPONENTIATES THE ACTION MATRIX (0,SIMAG)                             AAXV0878
C                         TO GIVE THE EDW S-MATRIX (SREAL,SIMAG)        AAXV0879
C NCHAN IS THE NUMBER OF CHANNELS                                       AAXV0880
C                                                                       AAXV0881
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0882
      DATA ZERO /0.0D0/                                                 AAXV0883
      DIMENSION SREAL(1),SIMAG(1),ETA(1),CA(NCHAN),SA(NCHAN),           AAXV0884
     ,                QMAT(NCHAN,NCHAN),SCR(NCHAN,5),ADIAG(NCHAN),      AAXV0885
     ,                                               VECTR(NCHAN,NCHAN) AAXV0886
C                                                                       AAXV0887
      DO 20 I  = 1,NCHAN                                                AAXV0888
      DO 10 J  = 1,5                                                    AAXV0889
  10  SCR(I,J) = ZERO                                                   AAXV0890
      DO 20 J  = 1,NCHAN                                                AAXV0891
  20  QMAT(I,J)= ZERO                                                   AAXV0892
C                                                                       AAXV0893
C DIAGONALIZE THE ACTION MATRIX                                         AAXV0894
C                                                                       AAXV0895
      K     = 0                                                         AAXV0896
      DO 30 I=1,NCHAN                                                   AAXV0897
      DO 30 J=1,I                                                       AAXV0898
      K     = K+1                                                       AAXV0899
      VECTR(I,J)=SIMAG(K)                                               AAXV0900
  30  CONTINUE                                                          AAXV0901
      CALL F02ABF(VECTR,NCHAN,NCHAN,ADIAG,QMAT,NCHAN,SCR,IFAIL)         AAXV0902
      DO 50 K = 1,NCHAN                                                 AAXV0903
      A = ADIAG(K)                                                      AAXV0904
      CA(K)   = DCOS(A)                                                 AAXV0905
  50  SA(K)   = DSIN(A)                                                 AAXV0906
C                                                                       AAXV0907
      N    = 0                                                          AAXV0908
      DO 300 N1 = 1,NCHAN                                               AAXV0909
      ETA1      = ETA(N1)                                               AAXV0910
      DO 200 N2 = 1,N1                                                  AAXV0911
      ETASUM    = ETA1+ETA(N2)                                          AAXV0912
      C  = DCOS(ETASUM)                                                 AAXV0913
      S  = DSIN(ETASUM)                                                 AAXV0914
      N  = N+1                                                          AAXV0915
      SR = ZERO                                                         AAXV0916
      SI = ZERO                                                         AAXV0917
C SUM OVER EIGENVALUES                                                  AAXV0918
      DO 100 K = 1,NCHAN                                                AAXV0919
      Q  = QMAT(N1,K)*QMAT(N2,K)                                        AAXV0920
      SR = SR+Q*CA(K)                                                   AAXV0921
      SI = SI+Q*SA(K)                                                   AAXV0922
  100 CONTINUE                                                          AAXV0923
      SREAL(N) = C*SR-S*SI                                              AAXV0924
      SIMAG(N) = S*SR+C*SI                                              AAXV0925
  200 CONTINUE                                                          AAXV0926
  300 CONTINUE                                                          AAXV0927
C                                                                       AAXV0928
      RETURN                                                            AAXV0929
      END                                                               AAXV0930
C                                                                       AAXV0931
      SUBROUTINE PROP(RN,RP,PPN,BBN,PPP,BBP,XK,ALFN,ALFP,QN,QP,WKBJ,    AAXV0932
     #                                                  CN,SN,CP,SP)    AAXV0933
C*******************************************************************    AAXV0934
C                                                                       AAXV0935
C PROPAGATES THE WAVEFUNCTION FROM RN TO RP AS DESCRIBED BY             AAXV0936
C       D.CHANG, J.C.LIGHT,,J.CHEM.PHYS. 50,2517 (1969)                 AAXV0937
C INPUT ARE THE SQUARES OF THE MOMENTA (TO BE CALLED POT.)              AAXV0938
C AND THEIR DERIVATIVES AT THE ENDPOINTS OF THE INTERVAL AS WELL AS     AAXV0939
C THE WAVEFUNCTION (S) AND DERIVATIVE (C) AND AMPLITUDE (ALF)           AAXV0940
C AND PHASE (Q) AT THE BEGINNING OF THE INTERVAL.                       AAXV0941
C OUTPUT ARE THESE LAST FOUR AT THE ENDPOINT.                           AAXV0942
C                                                                       AAXV0943
C                                                                       AAXV0944
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV0945
      LOGICAL WKBJ                                                      AAXV0946
C                                                                       AAXV0947
      DATA ZERO /0.0D0/                                                 AAXV0948
      DATA D2, D10 /0.5D0,0.1D0/                                        AAXV0949
      DATA  D6, D60/0.166666666666667D0,0.016666666666667D0/            AAXV0950
      DATA TWOPI /6.28318530717958D0/                                   AAXV0951
C                                                                       AAXV0952
      H = RP-RN                                                         AAXV0953
C                                                                       AAXV0954
      IF (WKBJ) THEN                                                    AAXV0955
C WKBJ - PROPAGATOR                                                     AAXV0956
C                                                                       AAXV0957
        WP = DSQRT(PPP)                                                 AAXV0958
        WN = DSQRT(PPN)                                                 AAXV0959
C                                                                       AAXV0960
C V = INTEGRAL OF SQRT(POT.) FROM RN TO RP                              AAXV0961
C                                                                       AAXV0962
        V  = H*(D2*(WP+WN) - H*(BBP/WP-BBN/WN)/24.0D0)                  AAXV0963
        S  = DSIN(V)                                                    AAXV0964
        C  = DCOS(V)                                                    AAXV0965
        SD = DSQRT(WP/WN)                                               AAXV0966
        SO = DSQRT(WP*WN)                                               AAXV0967
        GAA = C/SD                                                      AAXV0968
        GBB = C*SD                                                      AAXV0969
        GAB = S/SO                                                      AAXV0970
        GBA = -S*SO                                                     AAXV0971
        TT  = V                                                         AAXV0972
      ELSE                                                              AAXV0973
C SECOND MAGNUS APPROXIMATION                                           AAXV0974
C                                                                       AAXV0975
        DWN  = D2*H*BBN                                                 AAXV0976
        DWP  = D2*H*BBP                                                 AAXV0977
        WBAR = D2*(PPP+PPN)-D6*(DWP-DWN)                                AAXV0978
        WDER = D60*(DWP+DWN)-D10*(PPP-PPN)                              AAXV0979
        OBB  = H*H*WDER                                                 AAXV0980
        OAB  = H                                                        AAXV0981
        OBA  = -H*WBAR                                                  AAXV0982
        TT   = -OAB*OBA                                                 AAXV0983
        T    = TT-OBB*OBB                                               AAXV0984
        CALL RAY(T,C,S)                                                 AAXV0985
        GAB  = S*OAB                                                    AAXV0986
        GBA  = S*OBA                                                    AAXV0987
        GAA  = C-S*OBB                                                  AAXV0988
        GBB  = C+S*OBB                                                  AAXV0989
        IF (TT.LT.0.0D0)       TT = 0.0D0                               AAXV0990
        TT = DSQRT(TT)                                                  AAXV0991
      END IF                                                            AAXV0992
C                                                                       AAXV0993
C PROPAGATOR COEFFICIENTS HAVE NOW BEEN CALCULATED                      AAXV0994
C PROPAGATE WAVEFUNCTION (S) AND DERIVATIVE (C)                         AAXV0995
C                                                                       AAXV0996
      SP   = GAA*SN+GAB*CN                                              AAXV0997
      CP   = GBA*SN+GBB*CN                                              AAXV0998
C                                                                       AAXV0999
C COMPUTE NEW AMPLITUDE                                                 AAXV1000
C AMPLITUDE = SQRT(S**2+(C/K)**2)  (DEFINITION)                         AAXV1001
C                                                                       AAXV1002
      AFACT = DSQRT(SP**2+(CP/XK)**2)                                   AAXV1003
      ALFP  = ALFN*AFACT                                                AAXV1004
C                                                                       AAXV1005
C NORMALIZE SP/CP                                                       AAXV1006
C                                                                       AAXV1007
      SP   = SP/AFACT                                                   AAXV1008
      CP   = CP/AFACT                                                   AAXV1009
C                                                                       AAXV1010
C COMPUTE PHASE                                                         AAXV1011
C                                                                       AAXV1012
      QP = DATAN2(SP,CP/XK)                                             AAXV1013
      IF(QP.LT.ZERO)        QP = QP+TWOPI                               AAXV1014
      NG = IDINT((QN+TT)/TWOPI)                                         AAXV1015
      QP = QP+DFLOAT(NG)*TWOPI                                          AAXV1016
      IF (DABS(QN+TT-QP).GT.3.0D0) QP = QP + DSIGN(TWOPI,QN+TT-QP)      AAXV1017
C                                                                       AAXV1018
      RETURN                                                            AAXV1019
      END                                                               AAXV1020
C                                                                       AAXV1021
      REAL FUNCTION QUAD*8(JN,NTSTEP,JTSTEP)                            AAXV1022
C*******************************************************************    AAXV1023
C                                                                       AAXV1024
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1025
C                                                                       AAXV1026
C DETERMINE WEIGHT FOR THIS POINT ,USING VARIOUS QUADRATURE SCHEMES     AAXV1027
C SEE J.STOER, EINFUEHRUNG IN DIE NUMERISCHE MATHEMATIK I               AAXV1028
C        SPRINGER-VERLAG , BERLIN (1972) , P 102                        AAXV1029
C                                                                       AAXV1030
      JN    = JN + 1                                                    AAXV1031
      DJTOT = 1.0D0                                                     AAXV1032
      IF (NTSTEP.EQ.1.OR.JTSTEP.EQ.1)    GOTO 370                       AAXV1033
      DJTOT = DFLOAT(JTSTEP)                                            AAXV1034
      IF (MOD(NTSTEP,2).EQ.0.AND.(NTSTEP-JN).LT.4)      GOTO 350        AAXV1035
C                                                                       AAXV1036
C SIMPSON RULE                                                          AAXV1037
C                                                                       AAXV1038
      DJTOT = DJTOT/3.0D0                                               AAXV1039
      IF (JN.EQ.1.OR.JN.EQ.NTSTEP)       GOTO 370                       AAXV1040
      DJTOT = 2.0D0*DJTOT                                               AAXV1041
      IF (MOD(JN,2).NE.0)                GOTO 370                       AAXV1042
      DJTOT = 2.0D0*DJTOT                                               AAXV1043
      GOTO 370                                                          AAXV1044
C                                                                       AAXV1045
  350 IF (NTSTEP.NE.2) THEN                                             AAXV1046
        N = NTSTEP - JN + 1                                             AAXV1047
C                                                                       AAXV1048
C 3/8 RULE  ("PULCHERRIMA")                                             AAXV1049
C                                                                       AAXV1050
        IF (N.EQ.1) THEN                                                AAXV1051
          DJTOT = DJTOT*3.0D0/8.0D0                                     AAXV1052
C                                                                       AAXV1053
C DECENT 3/8 RULE                                                       AAXV1054
C                                                                       AAXV1055
        ELSE IF (N.EQ.2.OR.N.EQ.3) THEN                                 AAXV1056
          DJTOT = DJTOT*9.0D0/8.0D0                                     AAXV1057
C                                                                       AAXV1058
C 3/8 RULE ("PULCHERRIMA")                                              AAXV1059
C                                                                       AAXV1060
        ELSE IF (NTSTEP.EQ.4) THEN                                      AAXV1061
          DJTOT = DJTOT*3.0D0/8.0D0                                     AAXV1062
C                                                                       AAXV1063
C PARTLY SIMPSON AND PARTLY 3/8  : ==> 1/3 + 3/8                        AAXV1064
C                                                                       AAXV1065
        ELSE                                                            AAXV1066
          DJTOT = DJTOT * 17.0D0/24.0D0                                 AAXV1067
        END IF                                                          AAXV1068
      ELSE                                                              AAXV1069
C                                                                       AAXV1070
C TRAPEZIUM RULE                                                        AAXV1071
C                                                                       AAXV1072
        DJTOT = DJTOT/2.0D0                                             AAXV1073
      END IF                                                            AAXV1074
C                                                                       AAXV1075
  370 QUAD = DJTOT                                                      AAXV1076
C                                                                       AAXV1077
      RETURN                                                            AAXV1078
      END                                                               AAXV1079
C                                                                       AAXV1080
      SUBROUTINE START(RZMIN,RSTART,TP,WVEC,NROW,UCF,NBB,NVIB,MVIB,RR,  AAXV1081
     &                                                   POTRX,DPOTRX)  AAXV1082
C*********************************************************************  AAXV1083
C                                                                       AAXV1084
C INITIALISES THE FOLLOWING ARRAYS                                      AAXV1085
C TP(N1)       TURNING POINT1                                           AAXV1086
C RSTART(N1)   STARTING POINTS FOR INTEGRATION                          AAXV1087
C                                                                       AAXV1088
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1089
      DIMENSION RSTART(1),TP(1),WVEC(1),NROW(1),QPOS(1),NBB(1),UCF(1),  AAXV1090
     *    RR(1),POTRX(1),DPOTRX(1)                                      AAXV1091
C                                                                       AAXV1092
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV1093
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV1094
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV1095
     #              NVMIN,NVMAX                                         AAXV1096
      DATA ZERO/0.0D0/                                                  AAXV1097
      DATA RCH/0.001D0/                                                 AAXV1098
      DATA DELTA/1.0D-7/                                                AAXV1099
C                                                                       AAXV1100
      IF (IPRNT.GT.2)            WRITE(6,603)                           AAXV1101
C                                                                       AAXV1102
      RZMIN  = RMAX                                                     AAXV1103
      TP1    = RMIN                                                     AAXV1104
      FCH    = 21.2D0*DFACTOR                                           AAXV1105
      NC     = 0                                                        AAXV1106
      DO 130 N1 = 1,NCHAN                                               AAXV1107
      WVECN1 = WVEC(NROW(N1))                                           AAXV1108
      NC     = NC+NBB(N1)                                               AAXV1109
C FIND INNERMOST TURNING POINT                                          AAXV1110
      CALL TURN(RMIN,RMAX,DELTA,TP1,PP1,BB1,JFAIL,JCOUNT,MAXSTEP,       AAXV1111
     #         UCF(N1),WVECN1,NBB(N1),NVIB,MVIB,RR,POTRX,DPOTRX)        AAXV1112
      IF (BB1.LT.0.0D0)                 JFAIL=3                         AAXV1113
      IF(IPRNT .GT. 2 .OR. JFAIL .NE. 0)                                AAXV1114
     #                    WRITE(6,601) N1,JCOUNT,TP1,PP1,BB1,JFAIL      AAXV1115
      IF (JFAIL.NE.0) THEN                                              AAXV1116
        WRITE(6,609) MAXSTEP                                            AAXV1117
        STOP 11                                                         AAXV1118
      END IF                                                            AAXV1119
      TP(N1) = TP1                                                      AAXV1120
      RZZ    = TP1                                                      AAXV1121
      IF (TP1.GT.RMAX)                 GOTO 100                         AAXV1122
C SKIP THE REST IF WAVEFUNCTIONS ARE NOT REQUIRED                       AAXV1123
C RZ IS THE STARTING DISTANCE                                           AAXV1124
      RZ    = TP1                                                       AAXV1125
C FIND DISTANCE AT WHICH PHASE APPROX EQUALS QPMIN BY ITERATION         AAXV1126
C MOVE IN ONE WAVELENGTH FROM THE TURNING POINT                         AAXV1127
      WAVEL = 1.0D0/WVECN1                                              AAXV1128
      CHANGE= WAVEL                                                     AAXV1129
      RZZMIN= RZ-CHANGE                                                 AAXV1130
      RZZMAX= RZ                                                        AAXV1131
C                                                                       AAXV1132
C FIRST FIND R VALUE FOR WHICH (U-K**2)(TP1-R)**2 IS .GT.FCH            AAXV1133
C                                                                       AAXV1134
   50 CALL POTD(POTRX,DPOTRX,RR,W,DW,RZZMIN,UCF(N1),WVECN1,NBB(N1))     AAXV1135
      F = -W*(RZZMIN-RZ)**2                                             AAXV1136
      IF((F-FCH).LT.0.0D0) THEN                                         AAXV1137
        RZZMIN = RZZMIN-CHANGE                                          AAXV1138
        IF(RZZMIN.LT.DELTA)    RZZMIN=(RZZMIN+CHANGE)*0.5D0             AAXV1139
        GOTO 50                                                         AAXV1140
      END IF                                                            AAXV1141
   70 RZZ = (RZZMAX+RZZMIN)*0.5D0                                       AAXV1142
      CALL POTD(POTRX,DPOTRX,RR,W,DW,RZZ,UCF(N1),WVECN1,NBB(N1))        AAXV1143
      F  = -W*(RZZ-RZ)**2                                               AAXV1144
      IF((F-FCH).LT.0.0D0) THEN                                         AAXV1145
        RZZMAX = RZZ                                                    AAXV1146
        IF((RZZMAX-RZZMIN).LT.RCH)     GOTO 100                         AAXV1147
        GOTO 70                                                         AAXV1148
      END IF                                                            AAXV1149
      RZZMIN = RZZ                                                      AAXV1150
      IF((RZZMAX-RZZMIN).LT.RCH)       GOTO 100                         AAXV1151
      GOTO 70                                                           AAXV1152
  100 RZ  = RZZ                                                         AAXV1153
C                                                                       AAXV1154
C WE HAVE FOUND STARTING POINT                                          AAXV1155
      RSTART(N1) = RZ                                                   AAXV1156
      IF(RZ .LT. RZMIN)     RZMIN = RZ                                  AAXV1157
      IF(RZMIN .LE. ZERO .OR. RZMIN .GT. RMAX)    WRITE(6,600)RZ        AAXV1158
      IF (IPRNT.GT.2)        WRITE(6,602) RZ                            AAXV1159
      IF (TP(N1).GT.RMAX) THEN                                          AAXV1160
        RSTART(N1)= 2.0D0*RMAX                                          AAXV1161
        WRITE(6,610) N1                                                 AAXV1162
        TP(N1) = 0.0D0                                                  AAXV1163
      END IF                                                            AAXV1164
C                                                                       AAXV1165
  130 CONTINUE                                                          AAXV1166
      IF (RZMIN.LT.RMIN)    WRITE(6,611) RZMIN                          AAXV1167
C                                                                       AAXV1168
  600 FORMAT(' STARTING POINT FOR INTEGRATION= ',F10.5/)                AAXV1169
  601 FORMAT(1X,3H***,' TURNING POINT FOR CHANNEL ',I4,                 AAXV1170
     #       ' (FOUND IN ',I4,' STEPS)  : ',F15.8,                      AAXV1171
     #                 ' (W/DW ',2D12.4,'  /JFAIL ',I2,' )')            AAXV1172
  602 FORMAT(4X,'INTEGRATION STARTS AT : ',F15.8)                       AAXV1173
  603 FORMAT(///,15X,'-- START --',/)                                   AAXV1174
  609 FORMAT(1X,'WARNING: TURNING POINT NOT FOUND AFTER',I4,'STEPS')    AAXV1175
  610 FORMAT(1X,10(1H*),' TURNING POINT FOR CHANNEL ',I5,' GT RMAX ',   AAXV1176
     #                                ' ==> CHANNEL CLOSED ',10(1H*))   AAXV1177
  611 FORMAT(1X,10(1H$),' STARTING DISTANCE ( ',F15.8,' ) LT RMIN ',    AAXV1178
     #       10(1H$))                                                   AAXV1179
C                                                                       AAXV1180
      RETURN                                                            AAXV1181
      END                                                               AAXV1182
C                                                                       AAXV1183
      SUBROUTINE TURN(A,B,EPSIL,RMIN,WMIN,DWMIN,IFAIL,ICOUNT,IMAX,      AAXV1184
     #                        UCF,WV,NBB,NVIB,MVIB,RR,POTRX,DPOTRX)     AAXV1185
C*******************************************************************    AAXV1186
C                                                                       AAXV1187
C CALCULATES THE INNERMOST ZERO OF THE KINETIC ENERGY FUNCTION IN (A,B) AAXV1188
C TO WITHIN A DISTANCE EPSIL. FIRST STAGE IS TO SEARCH FOR A CHANGE     AAXV1189
C IN SIGN. SECOND STAGE IS TO CONTRACT THE INTERVAL CONTAINING          AAXV1190
C THE ZERO USING THE REGULA FALSI METHOD.                               AAXV1191
C        R=(WN*RP-WP*RN)/(WN-WP)                                        AAXV1192
C FINALLY USE NEWTON-RAPHSON TO OBTAIN THE PRECISE VALUE                AAXV1193
C THE ROUTINE IS SAFEGUARDED AGAINST A NOT WANTED TURNING POINT         AAXV1194
C IF THERE IS A REASONABLE LOW-LYING WELL LEFT OF THE                   AAXV1195
C TURNING POINT (I.E. IN THE FORBIDDEN REGION) THE                      AAXV1196
C ROUTINE CAN YIELD ITS BOTTOM AS TURNING POINT,                        AAXV1197
C SINCE THIS IS A PROPER STARTING POINT FOR SUBR. START                 AAXV1198
C                                                                       AAXV1199
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1200
      DIMENSION NBB(1),UCF(1),RR(1),POTRX(1),DPOTRX(1)                  AAXV1201
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV1202
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV1203
     #              NVMIN,NVMAX                                         AAXV1204
      DATA ZERO/0.0D0/                                                  AAXV1205
C                                                                       AAXV1206
C ESTIMATE REASONABLE STEPSIZE AND CLOSENESS CRITERIUM                  AAXV1207
      DELTA  = DMAX1((B-2.0D0*A)/50.0D0,A/3.0D0)                        AAXV1208
      EPSN   = DELTA/8.0D0                                              AAXV1209
C                                                                       AAXV1210
      IFAIL  = 0                                                        AAXV1211
      ICOUNT = 0                                                        AAXV1212
C TRY NEWTON METHOD                                                     AAXV1213
      CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)             AAXV1214
      IF(DABS(WMIN/DWMIN) .GE. DELTA) THEN                              AAXV1215
        RN   = A                                                        AAXV1216
        H    = DELTA                                                    AAXV1217
C LOOK FOR A CHANGE IN SIGN                                             AAXV1218
   50   CONTINUE                                                        AAXV1219
C STEP OUTWARDS                                                         AAXV1220
        ICOUNT= ICOUNT+1                                                AAXV1221
        RMIN  = RN+H                                                    AAXV1222
        CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)           AAXV1223
        IF(IPRNT .GT. 3)  WRITE(6,600)RMIN,WMIN,DWMIN,IFAIL,ICOUNT      AAXV1224
        IF(ICOUNT .GT. IMAX)         GOTO 300                           AAXV1225
        IF(WMIN .LT. ZERO)   THEN                                       AAXV1226
C NO CHANGE IN SIGN                                                     AAXV1227
          RN   = RMIN                                                   AAXV1228
          WN   = WMIN                                                   AAXV1229
          IF(RMIN .GT. B)            GOTO 300                           AAXV1230
          GOTO 50                                                       AAXV1231
        END IF                                                          AAXV1232
C KINETIC ENERGY CHANGES SIGN                                           AAXV1233
        RP   = RMIN                                                     AAXV1234
        WP   = WMIN                                                     AAXV1235
C REGULA FALSI METHOD                                                   AAXV1236
  150   CONTINUE                                                        AAXV1237
C FIND NEW MINIMUM POINT                                                AAXV1238
        Y   = (RP-RN)/(WP-WN)                                           AAXV1239
C FIND MINIMUM LENGTH OF TWO INTERVALS                                  AAXV1240
        WMIN=WP                                                         AAXV1241
        IF(-WN .LT. WP)       WMIN = -WN                                AAXV1242
        HMIN = WMIN*Y                                                   AAXV1243
        RMIN = RP-WP*Y                                                  AAXV1244
        CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)           AAXV1245
        IF(IPRNT .GT. 3)  WRITE(6,600)RMIN,WMIN,DWMIN,IFAIL,ICOUNT      AAXV1246
        ICOUNT = ICOUNT+1                                               AAXV1247
        IF(ICOUNT .GT. IMAX)             GOTO 300                       AAXV1248
        IF(HMIN .LT. EPSN)               GOTO 350                       AAXV1249
C STEP LENGTH STILL TOO LARGE                                           AAXV1250
        IF(WMIN .LE. ZERO)  THEN                                        AAXV1251
C NEGATIVE KINETIC ENERGY                                               AAXV1252
          RN   = RMIN                                                   AAXV1253
          WN   = WMIN                                                   AAXV1254
          GOTO 150                                                      AAXV1255
        END IF                                                          AAXV1256
C                                                                       AAXV1257
C POSITIVE KINETIC ENERGY                                               AAXV1258
        RP   = RMIN                                                     AAXV1259
        WP   = WMIN                                                     AAXV1260
        GOTO 150                                                        AAXV1261
  300   CONTINUE                                                        AAXV1262
C TOO MANY STEPS OR TURNING POINT TOO FAR OUT                           AAXV1263
        IF (RMIN.LE.B)       IFAIL=2                                    AAXV1264
        RETURN                                                          AAXV1265
      END IF                                                            AAXV1266
C                                                                       AAXV1267
C NEWTON METHOD                                                         AAXV1268
  350 IF (DWMIN.GE.0.0D0)  THEN                                         AAXV1269
        H = -WMIN/DWMIN                                                 AAXV1270
        IF (H.GT.DELTA)      H =  DELTA                                 AAXV1271
        IF (H.LT.-DELTA)     H = -DELTA                                 AAXV1272
C                                                                       AAXV1273
        RMIN = RMIN+H                                                   AAXV1274
        IF(RMIN.LT.0.01D0)        RMIN = 0.5D0*(RMIN-H)                 AAXV1275
        CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)           AAXV1276
        IF (IPRNT.GT.3)        WRITE(6,601) RMIN,WMIN,DWMIN,ICOUNT      AAXV1277
        ICOUNT = ICOUNT+1                                               AAXV1278
        IF(ICOUNT  .GT. IMAX)       GOTO 300                            AAXV1279
        IF(DABS(H) .GT. EPSIL)        GOTO 350                          AAXV1280
C NEWTON METHOD CONVERGED                                               AAXV1281
        RETURN                                                          AAXV1282
      END IF                                                            AAXV1283
C                                                                       AAXV1284
C DWMIN < 0   ==>  WRONG OR NON-EXISTING TURNING POINT                  AAXV1285
C                                                                       AAXV1286
      H    = DELTA                                                      AAXV1287
  920 RP   = RMIN                                                       AAXV1288
      RMIN = RMIN + H                                                   AAXV1289
      IF (ICOUNT.GT.IMAX .OR. RMIN.LT.0.0D0)     GOTO 300               AAXV1290
      CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)             AAXV1291
      IF (IPRNT.GT.3)       WRITE(6,602) RMIN,WMIN,DWMIN,ICOUNT         AAXV1292
      ICOUNT = ICOUNT + 1                                               AAXV1293
      IF (DWMIN.LT.0.0D0)            GOTO 920                           AAXV1294
C                                                                       AAXV1295
C DWMIN HAS CORRECT SIGN / CHECK WMIN                                   AAXV1296
C                                                                       AAXV1297
      IF (WMIN.LE.0.0D0)  THEN                                          AAXV1298
C                                                                       AAXV1299
C WMIN < 0  ==> GO TO BOTTOM                                            AAXV1300
C                                                                       AAXV1301
        IF (H.LT.EPSN)                GOTO 940                          AAXV1302
        H    = H/2.0D0                                                  AAXV1303
        RMIN = RP                                                       AAXV1304
        GOTO 920                                                        AAXV1305
      END IF                                                            AAXV1306
C                                                                       AAXV1307
C WMIN > 0  ==> GO BACK UNTIL SIGN CHANGES                              AAXV1308
C                                                                       AAXV1309
  930 RP   = RMIN                                                       AAXV1310
      WP   = WMIN                                                       AAXV1311
      RMIN = RMIN - H                                                   AAXV1312
      IF (RMIN.LT.0.0D0.OR.ICOUNT.GT.IMAX)      GOTO 300                AAXV1313
      CALL POTD(POTRX,DPOTRX,RR,WMIN,DWMIN,RMIN,UCF,WV,NBB)             AAXV1314
      IF (IPRNT.GT.3)       WRITE(6,603) RMIN,WMIN,DWMIN,ICOUNT         AAXV1315
      ICOUNT = ICOUNT + 1                                               AAXV1316
      RN   = RMIN                                                       AAXV1317
      WN   = WMIN                                                       AAXV1318
      IF (WN.GT.0.0D0)                GOTO 930                          AAXV1319
      GOTO 150                                                          AAXV1320
C                                                                       AAXV1321
  940 CONTINUE                                                          AAXV1322
      RETURN                                                            AAXV1323
C                                                                       AAXV1324
  600 FORMAT(' RMIN WMIN DWMIN IFAIL ICOUNT',3D20.12,2I6)               AAXV1325
  601 FORMAT(' (NR) : RMIN,WMIN,DWMIN,ICOUNT ',3D20.12,I6)              AAXV1326
  602 FORMAT(' (S1) : RMIN,WMIN,DWMIN,ICOUNT ',3D20.12,I6)              AAXV1327
  603 FORMAT(' (S2) : RMIN,WMIN,DWMIN,ICOUNT ',3D20.12,I6)              AAXV1328
C                                                                       AAXV1329
      END                                                               AAXV1330
      SUBROUTINE WRITEC(NJMIN1,NJMAX1,NJMI2,NJMA2,CONV,FLUX,XK,CRJ)     AAXV1331
C*******************************************************************    AAXV1332
C                                                                       AAXV1333
C PRINTS THE MATRIX OF CROSS SECTIONS                                   AAXV1334
C FOR INITIAL ENERGY LEVELS NJMIN1 TO NJMAX1                            AAXV1335
C AND FINAL ENERGY LEVELS NJMI2 TO NJMA2                                AAXV1336
C MULTIPLY WITH PROPER FACTORS IN ORDER TO GET THE DESIRED DIMENSIONED  AAXV1337
C (I.E. CONV*PI/(XK**2*(2J+1))) OR DIMENSIONLESS (I.E. 1/(2J+1))        AAXV1338
C CROSSECTION                                                           AAXV1339
C                                                                       AAXV1340
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1341
      DIMENSION CRJ(1),FLUX(1),XK(1),C(10)                              AAXV1342
C                                                                       AAXV1343
      DATA PI/3.14159265358979D0/                                       AAXV1344
C                                                                       AAXV1345
      KR  = (NJMA2 - NJMI2)/10 + 1                                      AAXV1346
      LRR = NJMA2 - NJMI2 + 1 - (KR-1)*10                               AAXV1347
      NJMAX2 = 0                                                        AAXV1348
      DO 300 K=1,KR                                                     AAXV1349
      LR = 10                                                           AAXV1350
      IF (K.EQ.KR) LR = LRR                                             AAXV1351
      NJMIN2 = NJMAX2 + 1                                               AAXV1352
      NJMAX2 = NJMAX2 + LR                                              AAXV1353
C                                                                       AAXV1354
      WRITE(6,602)     (NJ2,NJ2=NJMIN2,NJMAX2)                          AAXV1355
      WRITE(6,603)                                                      AAXV1356
      DO 200 NJ1 = NJMIN1,NJMAX1                                        AAXV1357
      FLUX1  = FLUX(NJ1)*CONV                                           AAXV1358
C                                                                       AAXV1359
C IF CONV.EQ.0 PRINT DIMENSIONLESS CROSSECTIONS                         AAXV1360
C                                                                       AAXV1361
      IF (DABS(CONV).LT.1.0D-20)       FLUX1 = FLUX(NJ1)*XK(NJ1)**2/PI  AAXV1362
      NJ = 0                                                            AAXV1363
      DO 100 NJ2=NJMIN2,NJMAX2                                          AAXV1364
      NJ = NJ + 1                                                       AAXV1365
      IF (NJ2.GT.NJ1) THEN                                              AAXV1366
        NJJ  = NJ2*(NJ2-1)/2+NJ1                                        AAXV1367
      ELSE                                                              AAXV1368
        NJJ  = NJ1*(NJ1-1)/2+NJ2                                        AAXV1369
      END IF                                                            AAXV1370
      C(NJ)  = FLUX1*CRJ(NJJ)                                           AAXV1371
  100 CONTINUE                                                          AAXV1372
  200 WRITE(6,600)    NJ1,(C(NJ2),NJ2=1,NJ)                             AAXV1373
  300 CONTINUE                                                          AAXV1374
      RETURN                                                            AAXV1375
C                                                                       AAXV1376
  600 FORMAT(1X,I4,5X,10D11.4)                                          AAXV1377
  602 FORMAT(1H0,/,7X,3HNV',10(I6,5X))                                  AAXV1378
  603 FORMAT(3X,2HNV)                                                   AAXV1379
C                                                                       AAXV1380
      END                                                               AAXV1381
C                                                                       AAXV1382
      SUBROUTINE ACTION(ETA,AMAT,XK,FK,NROW,UCF,NBB,                    AAXV1383
     #                  RSTART,SP,SN,CP,CN,ANEG,QNEG,APOS,QPOS,         AAXV1384
     #                        WN,WP,BB,TP,NVIB,MVIB,RR,POTRX,DPOTRX)    AAXV1385
C*******************************************************************    AAXV1386
C                                                                       AAXV1387
C CALCULATES A-MATRIX AND PHASE SHIFTS                                  AAXV1388
C                                                                       AAXV1389
C ETA(N1)             PHASE SHIFTS                                      AAXV1390
C AMAT(N)             ACTION MATRIX                                     AAXV1391
C APOS(N1),ANEG(N1)   AMPLITUDE OF CHANNEL WAVEFUNCTIONS                AAXV1392
C QPOS(N1),QNEG(N1)   PHASE OF WAVEFUNCTIONS                            AAXV1393
C SN,SP               CHANNEL WAVEFUNCTIONS                             AAXV1394
C CN,CP               DERIVATIVES OF CHANNEL WAVEFUNCTIONS              AAXV1395
C                                                                       AAXV1396
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1397
      LOGICAL WKBJ,SWKBJ                                                AAXV1398
C                                                                       AAXV1399
      COMMON/CHAN / RMIN,RMAX,DFACTOR,ACC,NCHAN,MAXSTEP,KJ              AAXV1400
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV1401
      COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ            AAXV1402
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV1403
     #              NVMIN,NVMAX                                         AAXV1404
C                                                                       AAXV1405
      DATA ZERO,HALF,POWER,HFMIN/0.0D0,0.5D0,0.25D0,0.2D0/              AAXV1406
      DATA APZ,APMAX/1.0D-20,1.0D30/                                    AAXV1407
C                                                                       AAXV1408
      DIMENSION ETA(1),AMAT(1),XK(1),FK(1),NROW(1),UCF(1),NBB(1),       AAXV1409
     #          RSTART(1),SP(1),SN(1),CP(1),CN(1),                      AAXV1410
     #          ANEG(1),QNEG(1),APOS(1),QPOS(1),WN(1),WP(1),            AAXV1411
     &          BB(1),TP(1),RR(1),POTRX(1),DPOTRX(1)                    AAXV1412
C                                                                       AAXV1413
C AUXILARY FUNCTION FOR DW-INTEGRALS (QQINT)                            AAXV1414
C                                                                       AAXV1415
      F(X,W,S,C) = (1.0D0-W/(X*X))*2.0D0*S*C                            AAXV1416
C                                                                       AAXV1417
      SWKBJ = .FALSE.                                                   AAXV1418
      WKBJ  = .FALSE.                                                   AAXV1419
C                                                                       AAXV1420
C FIND STARTING DISTANCES                                               AAXV1421
C                                                                       AAXV1422
      CALL START(RZMIN,RSTART,TP,XK,NROW,UCF,NBB,NVIB,MVIB,RR,          AAXV1423
     &                                            POTRX,DPOTRX)         AAXV1424
C                                                                       AAXV1425
      TPMAX    = 0.0D0                                                  AAXV1426
      DO 100 N1  = 1,NCHAN                                              AAXV1427
      APOS(N1) = 0.0D0                                                  AAXV1428
      QPOS(N1) = 0.0D0                                                  AAXV1429
      IF (TP(N1).GT.TPMAX)    TPMAX = TP(N1)                            AAXV1430
  100 CONTINUE                                                          AAXV1431
C                                                                       AAXV1432
      IF (IPRNT.GT.2)          WRITE(6,601)                             AAXV1433
      ICOUNT = 0                                                        AAXV1434
      H = ACC*100.0D0                                                   AAXV1435
      RPOS = RZMIN                                                      AAXV1436
      IF(IPRNT .GT. 3)         WRITE(6,604)                             AAXV1437
C                                                                       AAXV1438
C OPEN LOOP OVER DISTANCE                                               AAXV1439
C                                                                       AAXV1440
200   CONTINUE                                                          AAXV1441
      ICOUNT = ICOUNT+1                                                 AAXV1442
      NORM   = 0                                                        AAXV1443
      N  = 0                                                            AAXV1444
C INITIALISE ERRMAX TO THRESHOLD VALUE                                  AAXV1445
      ERRMAX = HFMIN*ACC                                                AAXV1446
      ERR    = ERRMAX                                                   AAXV1447
C                                                                       AAXV1448
C UPDATE ENDPOINTS AND CENTRE OF INTERVAL                               AAXV1449
C                                                                       AAXV1450
      RNEG  = RPOS                                                      AAXV1451
      RPOS  = RNEG+H                                                    AAXV1452
      RCEN  = RNEG+HALF*H                                               AAXV1453
C                                                                       AAXV1454
C INITIALIZE WKBJ - CHECK                                               AAXV1455
C                                                                       AAXV1456
      IF (.NOT.SWKBJ.AND.WKBJ.AND.IPRNT.GT.2)   WRITE(6,607) RNEG       AAXV1457
      IF (.NOT.WKBJ.AND.SWKBJ.AND.IPRNT.GT.2)   WRITE(6,606) RNEG       AAXV1458
      WKBJ  = .FALSE.                                                   AAXV1459
      IF (SWKBJ)    WKBJ = .TRUE.                                       AAXV1460
      SWKBJ = .TRUE.                                                    AAXV1461
      IF (LWKBJ.EQ.0 .OR. RNEG.LT.TPMAX)    SWKBJ = .FALSE.             AAXV1462
C                                                                       AAXV1463
C LOOP OVER INITIAL CHANNEL                                             AAXV1464
C                                                                       AAXV1465
      NC  = 0                                                           AAXV1466
      DO 260 N1=1,NCHAN                                                 AAXV1467
      RZ  = RSTART(N1)                                                  AAXV1468
      NC  = NC + NBB(N1)                                                AAXV1469
C SKIP IF INSIDE STARTING POINT                                         AAXV1470
      IF(RPOS .LT. RZ)                   GOTO 255                       AAXV1471
C                                                                       AAXV1472
      XK1 = XK(NROW(N1))                                                AAXV1473
C CHECK WHETHER STARTING POINT LIES WITHIN THIS INTERVAL                AAXV1474
      IF(RNEG .LE. RZ) THEN                                             AAXV1475
C                                                                       AAXV1476
C STARTING INTERVAL                                                     AAXV1477
C                                                                       AAXV1478
        CALL POTD(POTRX,DPOTRX,RR,WN1,BN1,RNEG,UCF(N1),XK1,NBB(N1))     AAXV1479
        CALL POTD(POTRX,DPOTRX,RR,WP1,BP1,RPOS,UCF(N1),XK1,NBB(N1))     AAXV1480
C                                                                       AAXV1481
C COMPUTE PHASE AND AMPLITUDE AT STARTING POINT RNEG                    AAXV1482
C                                                                       AAXV1483
        T = DSQRT(-WN1)-0.25D0*BN1/WN1                                  AAXV1484
        T = XK1/T                                                       AAXV1485
        QN1 = DATAN(T)                                                  AAXV1486
        AN1 = APZ                                                       AAXV1487
        IF (IPRNT.GT.2)       WRITE(6,609) N1,RNEG,QN1,AN1              AAXV1488
C                                                                       AAXV1489
C SET SINE AND COSINE OF PHASE AT STARTING POINT                        AAXV1490
C                                                                       AAXV1491
C SN IS THE WAVEFUNCTION                                                AAXV1492
C CN IS THE FIRST DERIVATIVE OF THE WAVEFUNCTION                        AAXV1493
C A IS THE AMPLITUDE  < = SQRT(SN**2+(CN/XK)**2)>                       AAXV1494
C Q IS THE PHASE      < = ATAN(SN/(CN/XK))>                             AAXV1495
C                                                                       AAXV1496
        CN(N1) = DCOS(QN1)*XK1                                          AAXV1497
        SN(N1) = DSIN(QN1)                                              AAXV1498
        CALL PROP(RNEG,RPOS,WN1,BN1,WP1,BP1,XK1,AN1,AE1,QN1,QE1,WKBJ,   AAXV1499
     #                                    CN(N1),SN(N1),CP(N1),SP(N1))  AAXV1500
        WN(N1) = WN1                                                    AAXV1501
        WP(N1) = WP1                                                    AAXV1502
        BB(N1) = BP1                                                    AAXV1503
        AP1    = AE1                                                    AAXV1504
        QP1    = QE1                                                    AAXV1505
        ERR    = ERRMAX                                                 AAXV1506
      ELSE                                                              AAXV1507
C                                                                       AAXV1508
C WAVEFUNCTIONS HAVE BEEN STARTED                                       AAXV1509
C                                                                       AAXV1510
C FIND WAVEFUNCTIONS FOR PREVIOUS INTERVAL                              AAXV1511
C                                                                       AAXV1512
        AN1 = APOS(N1)                                                  AAXV1513
        QN1 = QPOS(N1)                                                  AAXV1514
C                                                                       AAXV1515
        WN1 = WP(N1)                                                    AAXV1516
        BN1 = BB(N1)                                                    AAXV1517
        CALL POTD(POTRX,DPOTRX,RR,WP1,BP1,RPOS,UCF(N1),XK1,NBB(N1))     AAXV1518
        WN(N1) = WN1                                                    AAXV1519
        WP(N1) = WP1                                                    AAXV1520
        BB(N1) = BP1                                                    AAXV1521
C                                                                       AAXV1522
C CHECK ON WKBJ CONDITION ACCORDING TO                                  AAXV1523
C D.CHANG, J.C.LIGHT,  J. CHEM. PHYS. 50, 2517 (1969)                   AAXV1524
C                                                                       AAXV1525
        IF (SWKBJ) THEN                                                 AAXV1526
          IF (WP1.LT.1.0D-30) THEN                                      AAXV1527
            SQWP1 = DSQRT(WP1)                                          AAXV1528
            IF (BP1.LT.BN1)        GOTO 226                             AAXV1529
            IF (TPMAX*SQWP1.LT.10.0D0)      GOTO 227                    AAXV1530
  226       IF (DABS(BP1)/(4.0D0*WP1*SQWP1).LT.STWKBJ)   GOTO 228       AAXV1531
          END IF                                                        AAXV1532
  227     SWKBJ = .FALSE.                                               AAXV1533
        END IF                                                          AAXV1534
C                                                                       AAXV1535
C SET SINE AND COSINE AT START OF PRESENT INTERVAL EQUAL TO THOSE AT    AAXV1536
C THE END OF PREVIOUS INTERVAL                                          AAXV1537
C                                                                       AAXV1538
  228   CN(N1) = CP(N1)                                                 AAXV1539
        SN(N1) = SP(N1)                                                 AAXV1540
C                                                                       AAXV1541
C WAVEFUNCTIONS PROPAGATE OVER WHOLE INTERVAL                           AAXV1542
C                                                                       AAXV1543
        CALL PROP(RNEG,RPOS,WN1,BN1,WP1,BP1,XK1,AN1,AE1,QN1,QE1,WKBJ,   AAXV1544
     #                                   CN(N1),SN(N1),CP(N1),SP(N1))   AAXV1545
C                                                                       AAXV1546
C NOW WAVEFUNCTIONS PROPAGATE OVER TWO HALF INTERVALS                   AAXV1547
C                                                                       AAXV1548
        CALL POTD(POTRX,DPOTRX,RR,WC1,BC1,RCEN,UCF(N1),XK1,NBB(N1))     AAXV1549
C                                                                       AAXV1550
        CALL PROP(RNEG,RCEN,WN1,BN1,WC1,BC1,XK1,AN1,AC1,QN1,QC1,WKBJ,   AAXV1551
     #                                       CN(N1),SN(N1),SVCN,SVSN)   AAXV1552
C                                                                       AAXV1553
        CALL PROP(RCEN,RPOS,WC1,BC1,WP1,BP1,XK1,AC1,AP1,QC1,QP1,WKBJ,   AAXV1554
     #                                       SVCN,SVSN,CP(N1),SP(N1))   AAXV1555
C                                                                       AAXV1556
C ESTIMATE ERROR FROM DIFFERENCE OF TWO PHASES                          AAXV1557
C                                                                       AAXV1558
        ERR = DABS(QP1-QE1)                                             AAXV1559
        AE1 = AP1                                                       AAXV1560
        QE1 = QP1                                                       AAXV1561
      END IF                                                            AAXV1562
C                                                                       AAXV1563
C WAVEFUNCTIONS STORED AT END POINTS OF INTERVAL                        AAXV1564
C (USE RESULTS FROM TWO HALF INTERVALS)                                 AAXV1565
C                                                                       AAXV1566
      ANEG(N1) = AN1                                                    AAXV1567
      QNEG(N1) = QN1                                                    AAXV1568
      APOS(N1) = AE1                                                    AAXV1569
      QPOS(N1) = QE1                                                    AAXV1570
      IF(AE1.GT.APMAX)          NORM=1                                  AAXV1571
      IF (IPRNT.GE.4) THEN                                              AAXV1572
C                                                                       AAXV1573
C PRINT RESULTS FOR THIS PROPAGATION STEP                               AAXV1574
C ( THE PHASE-SHIFT WILL NOT BE VERY ACCURATE)                          AAXV1575
C                                                                       AAXV1576
C FOR PHASE-CALC. QRB IS CALLED WITH RL CONSISTENT WITH                 AAXV1577
C RL*(RL+1) = UCF(N1)                                                   AAXV1578
C                                                                       AAXV1579
        RL = (DSQRT(1.0D0+4.0D0*UCF(N1)) - 1.0D0)/2.0D0                 AAXV1580
C                                                                       AAXV1581
        ETA1 = QP1-QRB(RL,XK1*RPOS)                                     AAXV1582
        DA   = (AP1-AN1)/(H*AP1)                                        AAXV1583
        DQ   = (QP1-QN1)/H-XK1                                          AAXV1584
        WRITE(6,602) ICOUNT,N1,RPOS,RNEG,AP1,QP1,DA,DQ,WP1,ERR,ETA1     AAXV1585
      END IF                                                            AAXV1586
C                                                                       AAXV1587
C END OF PROPAGATION STEP                                               AAXV1588
C                                                                       AAXV1589
      IF(ERR .GT. ERRMAX)      ERRMAX = ERR                             AAXV1590
      IF(IDWINTG .GT. 0) THEN                                           AAXV1591
C LOOP OVER FINAL CHANNEL                                               AAXV1592
        NB = N1  - NBB(N1) + 1                                          AAXV1593
        N  = N + NB                                                     AAXV1594
        NE = N1 - 1                                                     AAXV1595
        IF (NE.LT.NB)                     GOTO 260                      AAXV1596
C       N1 A-MATRIX CONTRIBUTIONS                                       AAXV1597
        FN1 = F(XK1,WN1,SN(N1),CN(N1))                                  AAXV1598
        FP1 = F(XK1,WP1,SP(N1),CP(N1))                                  AAXV1599
        DO 250 N2 = NB,NE                                               AAXV1600
C SKIP IF IN FORBIDDEN REGION                                           AAXV1601
        IF (RPOS.GE.RSTART(N2)) THEN                                    AAXV1602
          K = N2 - N1                                                   AAXV1603
C                                                                       AAXV1604
C SKIP IF ZERO-CONTRIBUTION                                             AAXV1605
          CALL POTND(POTRX,DPOTRX,RR,WP21,BP21,RPOS,UCF(N1),XK1,N1,N2)  AAXV1606
          CALL POTND(POTRX,DPOTRX,RR,WN21,BN21,RNEG,UCF(N1),XK1,N1,N2)  AAXV1607
C                                                                       AAXV1608
C CONTRIBUTION TO ACTION MATRIX OVER THIS INTERVAL                      AAXV1609
C                                                                       AAXV1610
          FN2 = F(XK(NROW(N2)),WN(N2),SN(N2),CN(N2))                    AAXV1611
          FP2 = F(XK(NROW(N2)),WP(N2),SP(N2),CP(N2))                    AAXV1612
          QQ  = QQINT(RNEG,RPOS,                                        AAXV1613
     #                ANEG(N2),APOS(N2),QNEG(N2),QPOS(N2),FN2,FP2,      AAXV1614
     #                    AN1,AP1,QN1,QP1,FN1,FP1,WN21,WP21,BN21,BP21)  AAXV1615
C UPDATE ACTION MATRIX                                                  AAXV1616
          AMAT(N) = AMAT(N)-QQ                                          AAXV1617
C                                                                       AAXV1618
          IF(IPRNT .GE. 4) THEN                                         AAXV1619
C PRINT WAVE FUNCTIONS AND ACTION MATRIX DENSITY                        AAXV1620
            QH = QQ/H                                                   AAXV1621
            WRITE(6,603)      N1,N2,WP21,QQ,QH,AMAT(N)                  AAXV1622
C CLOSE LOOP ON FINAL CHANNEL                                           AAXV1623
          END IF                                                        AAXV1624
        END IF                                                          AAXV1625
  250   N = N + 1                                                       AAXV1626
        GOTO 260                                                        AAXV1627
C ADJUST N (FOR RPOS.LT.RSTART)                                         AAXV1628
      END IF                                                            AAXV1629
  255 N = N + N1                                                        AAXV1630
C                                                                       AAXV1631
C CLOSE LOOP ON INITIAL CHANNEL                                         AAXV1632
C                                                                       AAXV1633
  260 CONTINUE                                                          AAXV1634
C                                                                       AAXV1635
      IF(NORM.GT.0)       CALL RENORM(NCHAN,APOS,ANEG,AMAT)             AAXV1636
C                                                                       AAXV1637
C CHANGE STEP LENGTH                                                    AAXV1638
C                                                                       AAXV1639
      HFAC  = (ACC/ERRMAX)**POWER                                       AAXV1640
      IF(HFAC .LT. HFMIN)    HFAC = HFMIN                               AAXV1641
      H  = HFAC*H                                                       AAXV1642
C                                                                       AAXV1643
C CHECK WHETHER STILL INSIDE INTERACTION REGION                         AAXV1644
C                                                                       AAXV1645
      IF(RPOS.LT.RMAX-H)                GOTO 200                        AAXV1646
C                                                                       AAXV1647
C CALCULATE PHASE SHIFTS AND NORMALISE ACTION MATRIX                    AAXV1648
C                                                                       AAXV1649
      N  = 0                                                            AAXV1650
      DO 310 N1=1,NCHAN                                                 AAXV1651
      XK1 = XK(NROW(N1))                                                AAXV1652
      QP1 = QPOS(N1)                                                    AAXV1653
      AP1 = APOS(N1)                                                    AAXV1654
C                                                                       AAXV1655
C CALCULATE PHASE SHIFTS                                                AAXV1656
C                                                                       AAXV1657
C FOR PHASE-CALC. QRB IS CALLED WITH RL CONSISTENT WITH                 AAXV1658
C RL*(RL+1) = UCF(N1)                                                   AAXV1659
C                                                                       AAXV1660
      RL = (DSQRT(1.0D0+4.0D0*UCF(N1)) - 1.0D0)/2.0D0                   AAXV1661
      ETA(N1) = QP1-QRB(RL,XK1*RPOS)                                    AAXV1662
      IF (RSTART(N1).GT.RMAX)   ETA(N1) = 0.0D0                         AAXV1663
      IF (IPRNT.GT.2)            WRITE(6,610) N1,RPOS,QP1,AP1           AAXV1664
C                                                                       AAXV1665
C NORMALIZE ACTION MATRIX                                               AAXV1666
C                                                                       AAXV1667
      T1 = ZERO                                                         AAXV1668
      IF(AP1 .GT. ZERO)       T1 = FK(NROW(N1))/AP1                     AAXV1669
      NB = N1 - NBB(N1)  + 1                                            AAXV1670
      NE = N1 - 1                                                       AAXV1671
      N  = N + NB                                                       AAXV1672
      IF (NE.GE.NB) THEN                                                AAXV1673
        DO 300 N2 = NB,NE                                               AAXV1674
        AP2 = APOS(N2)                                                  AAXV1675
        IF (AP2.GT.ZERO)   AMAT(N)=AMAT(N)*T1*FK(NROW(N2))/AP2          AAXV1676
  300   N = N + 1                                                       AAXV1677
      END IF                                                            AAXV1678
  310 CONTINUE                                                          AAXV1679
C                                                                       AAXV1680
      IF(IPRNT.GT.1)        WRITE(6,605) RZMIN,RPOS,ICOUNT              AAXV1681
C                                                                       AAXV1682
  601 FORMAT(///,15X,'-- ACTION --',/,1X)                               AAXV1683
  605 FORMAT(1H0,//,40(1H.),' ACTION : ',                               AAXV1684
     #'A-MATRIX FROM',F10.5,' TO',F10.5,' IN',I7,' STEPS',/)            AAXV1685
  602 FORMAT(' W',2I4,F10.5,7D12.4,F10.5)                               AAXV1686
  603 FORMAT(' C',2I4,3D12.4,D15.7)                                     AAXV1687
  604 FORMAT(//,1X,2H**,' W PRINTS : NUMBER/CHANNEL/  R /STEPSIZE/'     AAXV1688
     #       ,'AMPLITUDE/PHASE/DEL AMPL/DEL PHAS/POT.SLOPE/ERROR/'      AAXV1689
     #       ,'PHASESHIFT',/,14X,  ' ICOUNT    N1   RPOS      H     '   AAXV1690
     #       ,'AP1     QP1     DA       DQ      PPP       ERR '         AAXV1691
     #       ,'   ETA1',//,1X,2H**,' C        : CHANNEL1/CHANNEL2'      AAXV1692
     #       ,'/INT.POT/A-MATRIX-CONTR./NORMALIZED A-MAT-C./'           AAXV1693
     #       ,'A-MATRIX-ELEMENT',//,1X)                                 AAXV1694
  609 FORMAT(4X,' CHANNEL ',I4,' STARTS AT R = ',F15.8,                 AAXV1695
     #       ' WITH PHASE ',D14.6,' AND AMPLITUDE ',D14.6)              AAXV1696
  610 FORMAT(4X,' CHANNEL ',I4,'  STOPS AT R = ',F15.8,                 AAXV1697
     #       ' WITH PHASE ',D14.6,' AND AMPLITUDE ',D14.6)              AAXV1698
  606 FORMAT(1H0,4X,5(1H*),' WKBJ PROPAGATION STARTS AT : ',F15.8,/,1X) AAXV1699
  607 FORMAT(1H0,4X,5(1H/),' WKBJ PROPAGATION STOPS AT :  ',F15.8,/,1X) AAXV1700
C                                                                       AAXV1701
      RETURN                                                            AAXV1702
      END                                                               AAXV1703
C                                                                       AAXV1704
      REAL FUNCTION QQINT*8(RN,RP,AN2,AP2,QN2,QP2,FN2,FP2,              AAXV1705
     #                      AN1,AP1,QN1,QP1,FN1,FP1,VN,VP,DVN,DVP)      AAXV1706
C*******************************************************************    AAXV1707
C                                                                       AAXV1708
C -- COMPUTATION OF DISTORTED WAVE INTEGRAL --                          AAXV1709
C INTEGRAL FROM RN TO RP OF                                             AAXV1710
C 2*A2(R)*V(R)*A1(R)*SIN(Q2(R))*SIN(Q1(R))                              AAXV1711
C THESE FUNCTIONS ARE INTERPOLATING POLYNOMIALS WITH VALUES             AAXV1712
C AJ(R)   ANJ,APJ                                                       AAXV1713
C QJ(R)   QNJ,QPJ                                                       AAXV1714
C WJ(R)   WNJ,WPJ                                                       AAXV1715
C V(R)    VN,VP,DVN,DVP                                                 AAXV1716
C AJ(R)*SIN(QJ(R)) ARE SOLUTIONS OF THE DIFFERENTIAL EQUATION           AAXV1717
C (D2/DR2 +WJ(R))AJ(R)*SIN(QJ(R))=0                                     AAXV1718
C WHERE WJ(R) HAS THE ASYMPTOTIC VALUE XKJ**2                           AAXV1719
C                                                                       AAXV1720
C F = (1-W/X**2)*2*SIN*COS IS PRECOMPUTED IN THE CALLING ROUTINE        AAXV1721
C (ACTION) AND GIVEN AS FN AND FP                                       AAXV1722
C                                                                       AAXV1723
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1724
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV1725
C                                                                       AAXV1726
      H4   = 0.25D0*(RP-RN)                                             AAXV1727
      AP   = AP2*AP1                                                    AAXV1728
      AN   = AN2*AN1                                                    AAXV1729
      GPOS = VP*AP+VN*AN                                                AAXV1730
      GNEG = VP*AP-VN*AN                                                AAXV1731
      DGN  = H4*AN*(DVN+DVN+(FN2+FN1)*VN)                               AAXV1732
      DGP  = H4*AP*(DVP+DVP+(FP2+FP1)*VP)                               AAXV1733
      DGPOS = DGP+DGN                                                   AAXV1734
      DGNEG = DGP-DGN                                                   AAXV1735
      PPOS = 0.5D0*(QP2-QN2+QP1-QN1)                                    AAXV1736
      PNEG = 0.5D0*(QP2-QN2-QP1+QN1)                                    AAXV1737
      QPOS = 0.5D0*(QP2+QN2+QP1+QN1)                                    AAXV1738
      QNEG = 0.5D0*(QP2+QN2-QP1-QN1)                                    AAXV1739
      CALL CAN(PPOS,C0,C1,C2,C3)                                        AAXV1740
      QSUM = DCOS(QPOS)*(GPOS*(C0+C0)+DGNEG*(C2-C0))                    AAXV1741
     +           +DSIN(QPOS)*(GNEG*(C3-3.0D0*C1)+DGPOS*(C1-C3))         AAXV1742
      CALL CAN(PNEG,C0,C1,C2,C3)                                        AAXV1743
      QDIF = DCOS(QNEG)*(GPOS*(C0+C0)+DGNEG*(C2-C0))                    AAXV1744
     +           +DSIN(QNEG)*(GNEG*(C3-3.0D0*C1)+DGPOS*(C1-C3))         AAXV1745
      QQINT = H4*(QSUM-QDIF)                                            AAXV1746
C                                                                       AAXV1747
      IF(IPRNT .LT. 6)       RETURN                                     AAXV1748
      WRITE(6,600) AP,PPOS,PNEG,VP,QDIF,QSUM                            AAXV1749
  600 FORMAT(1H ,6D15.7)                                                AAXV1750
      RETURN                                                            AAXV1751
      END                                                               AAXV1752
C                                                                       AAXV1753
      SUBROUTINE RAY(T,C,S)                                             AAXV1754
C*******************************************************************    AAXV1755
C                                                                       AAXV1756
C  C(T)  = COS(SQRT(T))                                                 AAXV1757
C  S(T)  = SIN(SQRT(T))/SQRT(T)                                         AAXV1758
C                                                                       AAXV1759
C CALLED BY PROP                                                        AAXV1760
C                                                                       AAXV1761
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1762
      IF(T .LT. -2500D0) THEN                                           AAXV1763
        WRITE(6,600)  T                                                 AAXV1764
        T  = -2500D0                                                    AAXV1765
      END IF                                                            AAXV1766
      IF(T.LE.-1.0D-3) THEN                                             AAXV1767
C T LESS THAN -0.001                                                    AAXV1768
        RT = DSQRT(-T)                                                  AAXV1769
        E  = DEXP(RT)                                                   AAXV1770
        C  = 0.5D0*E+0.5D0/E                                            AAXV1771
        S  = (0.5D0*E-0.5D0/E)/RT                                       AAXV1772
      ELSEIF(T.LE.1.0D-3) THEN                                          AAXV1773
C T VERY SMALL                                                          AAXV1774
        C = 1.0D0+T*(T-12.0D0)/24.0D0                                   AAXV1775
        S = 1.0D0+T*(T-20.0D0)/120.0D0                                  AAXV1776
      ELSE                                                              AAXV1777
C T GREATER THAN 0.001                                                  AAXV1778
        RT = DSQRT(T)                                                   AAXV1779
        C  = DCOS(RT)                                                   AAXV1780
        S  = DSIN(RT)/RT                                                AAXV1781
      END IF                                                            AAXV1782
      RETURN                                                            AAXV1783
C                                                                       AAXV1784
  600 FORMAT(' WARNING: ARGUMENT OF RAY TOO LARGE T=',D20.10)           AAXV1785
      END                                                               AAXV1786
C                                                                       AAXV1787
      INTEGER FUNCTION MAINCS(AR,MAXC,NBYTE)                            AAXV1788
C*******************************************************************    AAXV1789
C                                                                       AAXV1790
C ROUTINE TO IMITATE DYNAMICAL CORE-REQUEST                             AAXV1791
C THE ARRAY (AR) REQUESTED IS SUPPOSED TO BE IN BLANK COMMON            AAXV1792
C STARTING AT THE FIRST ADDRESS THEREIN                                 AAXV1793
C (NOTE:...THE ARRAY AR IS NOT USED OR CHECKED...)                      AAXV1794
C THIS ROUTINE MERELY CHECKS IF NCORE IS LARGE ENOUGH                   AAXV1795
C (NCORE IS SET EQUAL TO THE DIMENSION OF BLANK COMMON                  AAXV1796
C IN THE MAIN PROGRAM)                                                  AAXV1797
C NOTE: NCORE SHOULD BE AT LEAST 1 MORE THEN THE NUMBER                 AAXV1798
C OF WORDS NEEDED                                                       AAXV1799
C                                                                       AAXV1800
C THE FOLLOWING VALUES ARE RETURNED:                                    AAXV1801
C ENOUGH CORE AVAILABLE          :   1                                  AAXV1802
C NOT ENOUGH CORE                :   0                                  AAXV1803
C NBYTE NOT EQUAL TO 8           :  -5                                  AAXV1804
C                                                                       AAXV1805
      IMPLICIT INTEGER (A-Z)                                            AAXV1806
C                                                                       AAXV1807
      COMMON/CORE / KSCR,KSR,KSI,KETA,KCRJ,KTOTAL,KXK,KFK,KCRJG,        AAXV1808
     $             NCHM,NCORE,KCR,KRR,KPOTRX,KDPOTRX                    AAXV1809
      K     = MAXC + 1                                                  AAXV1810
      IF (NBYTE.NE.8) THEN                                              AAXV1811
        MAINCS = -5                                                     AAXV1812
        RETURN                                                          AAXV1813
      ELSE                                                              AAXV1814
        MAINCS =  1                                                     AAXV1815
      END IF                                                            AAXV1816
      IF(K.GT.NCORE) THEN                                               AAXV1817
        MAINCS = 0                                                      AAXV1818
        WRITE(6,601) K                                                  AAXV1819
      END IF                                                            AAXV1820
      RETURN                                                            AAXV1821
C                                                                       AAXV1822
  601 FORMAT(1H0,10(1H*),'NOT ENOUGH BLANK COMMON AVAILABLE',/,         AAXV1823
     $       1X,10(1H*),' NCORE AND THE DIMENSION OF CR IN THE MAIN'    AAXV1824
     $       ,' PROGRAM SHOULD BE AT LEAST : ',I6,/,1H0)                AAXV1825
C                                                                       AAXV1826
      END                                                               AAXV1827
C                                                                       AAXV1828
      SUBROUTINE WRITER(A,NTOT)                                         AAXV1829
C*******************************************************************    AAXV1830
C                                                                       AAXV1831
C PRINT A TRIANGULAR MATRIX                                             AAXV1832
C                                                                       AAXV1833
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1834
      DIMENSION A(1)                                                    AAXV1835
C                                                                       AAXV1836
      DO 101 L=1,NTOT,10                                                AAXV1837
      WRITE(6,108)                                                      AAXV1838
      IEND = MIN0(NTOT,L+10-1)                                          AAXV1839
      WRITE(6,107)   (I,I=L,IEND)                                       AAXV1840
      WRITE(6,608)                                                      AAXV1841
      DO 103 K=L,NTOT                                                   AAXV1842
      KL = (K*(K-1))/2+L                                                AAXV1843
      KM = KL + MIN0(K-L,9)                                             AAXV1844
 103  WRITE(6,104)    K,(A(I),I=KL,KM)                                  AAXV1845
 101  CONTINUE                                                          AAXV1846
      RETURN                                                            AAXV1847
C                                                                       AAXV1848
 104  FORMAT(1X,I3,2X,10D11.4)                                          AAXV1849
 107  FORMAT(5X,I6,9I11)                                                AAXV1850
 108  FORMAT(/)                                                         AAXV1851
 608  FORMAT(1X)                                                        AAXV1852
C                                                                       AAXV1853
      END                                                               AAXV1854
C                                                                       AAXV1855
      SUBROUTINE POTMAT(RR,GAMMA,XVALUE,POTRX,DPOTRX,NPOINT,VIBM)       AAXV1856
C*******************************************************************    AAXV1857
C EVALUATES POTENTIAL MATRIX, GIVEN WAVE FUNCTIONS AND POTENTIAL        AAXV1858
C                                                                       AAXV1859
      IMPLICIT REAL*8 (A-H,O-Z)                                         AAXV1860
      PARAMETER(NDIM=96,N2DIM=4)                                        AAXV1861
      COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEP,   AAXV1862
     #              NVMIN,NVMAX                                         AAXV1863
      COMMON/OPT  / IPRNT,IDWINTG,IDISC                                 AAXV1864
      DIMENSION XVALUE(NDIM),VIBM(NDIM,N2DIM),FUNC(NDIM),DFUNC(NDIM),   AAXV1865
     *         VIBPHI(NDIM),VIBLAG(NDIM),RR(1),POTRX(1),DPOTRX(1)       AAXV1866
C                                                                       AAXV1867
          IFAIL = 0                                                     AAXV1868
          III   = 1                                                     AAXV1869
          IPOINT = NDIM                                                 AAXV1870
          NVMXP1 = N2DIM                                                AAXV1871
C                                                                       AAXV1872
          DO 300 IR=1,NPOINT                                            AAXV1873
          DO 200 NVIB=1,NVMXP1                                          AAXV1874
          DO 160 MVIB=1,NVIB                                            AAXV1875
          DO 125 LX=1,IPOINT                                            AAXV1876
          RI = XVALUE(LX)                                               AAXV1877
          CALL POTENT(POT,DPOT,GAMMA,RR(IR),RI)                         AAXV1878
          FUNC(LX)  = VIBM(LX,MVIB)*VIBM(LX,NVIB)*POT                   AAXV1879
          DFUNC(LX) = VIBM(LX,MVIB)*VIBM(LX,NVIB)*DPOT                  AAXV1880
  125     CONTINUE                                                      AAXV1881
          CALL D01GAF(XVALUE,FUNC,IPOINT,VINTG,ERR,IFAIL)               AAXV1882
          CALL D01GAF(XVALUE,DFUNC,IPOINT,DVINTG,ERR,IFAIL)             AAXV1883
          POTRX(III)  = VINTG                                           AAXV1884
          DPOTRX(III) = DVINTG                                          AAXV1885
          III = III+1                                                   AAXV1886
  160     CONTINUE                                                      AAXV1887
  200     CONTINUE                                                      AAXV1888
  300     CONTINUE                                                      AAXV1889
          RETURN                                                        AAXV1890
          END                                                           AAXV1891
C                                                                       AAXV1892
        SUBROUTINE POTENT(POT,DPOT,GAMMA,R,RI)                          AAXV1893
C*******************************************************************    AAXV1894
C EVALUATES POTENTIAL USING THE MORSE FUNCTION                          AAXV1895
C                                                                       AAXV1896
         IMPLICIT REAL*8(A-H,O-Z)                                       AAXV1897
         COMMON/PARAM/ RVDW,DVDW,ALVDW,RDIA,DDIA,ALDIA,ZXK              AAXV1898
C                                                                       AAXV1899
         R0   = RVDW                                                    AAXV1900
         D    = DVDW                                                    AAXV1901
         ALPHA = ALVDW                                                  AAXV1902
         RB   = DSQRT(R*R+RI*RI/4.-RI*R*DCOS(GAMMA))                    AAXV1903
         RA   = DSQRT(R*R+RI*RI/4.+RI*R*DCOS(GAMMA))                    AAXV1904
         EXRA = DEXP(-ALPHA*(RA-R0))                                    AAXV1905
         EXRB = DEXP(-ALPHA*(RB-R0))                                    AAXV1906
         POT  = D*((EXRA-1.D0)*(EXRA-1.D0)+(EXRB-1.D0)*(EXRB-1.D0)-2.D0)AAXV1907
         DPOT =-ALPHA*D*((EXRA-1.D0)*EXRA*(2.D0*R+RI*                   AAXV1908
     *               DCOS(GAMMA))/RA+(EXRB-1.D0)*EXRB*(2.D0*R-          AAXV1909
     -                                          RI*DCOS(GAMMA))/RB)     AAXV1910
         RETURN                                                         AAXV1911
         END                                                            AAXV1912
C                                                                       AAXV1913
        SUBROUTINE VIBR(XVALUE,VIBM)                                    AAXV1914
C*******************************************************************    AAXV1915
C CALCULATES THE VIBRATIONAL WAVE FUNCTIONS                             AAXV1916
C                                                                       AAXV1917
          IMPLICIT REAL*8(A-H,O-Z)                                      AAXV1918
          REAL RINDEX                                                   AAXV1919
          PARAMETER(NDIM=96,N2DIM=4)                                    AAXV1920
          COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ        AAXV1921
          COMMON/OPT  / IPRNT,IDWINTG,IDISC                             AAXV1922
          COMMON/PARAM/ RVDW,DVDW,ALVDW,RDIA,DDIA,ALDIA,ZXK             AAXV1923
          DIMENSION VIBLAG(NDIM),VIBPHI(NDIM),XVALUE(NDIM),             AAXV1924
     &                                          VIBM(NDIM,N2DIM)        AAXV1925
C                                                                       AAXV1926
          IFAIL  = 0                                                    AAXV1927
          NP     = NDIM                                                 AAXV1928
          NVMXP1 = N2DIM                                                AAXV1929
          VINDEX = 0.0D0                                                AAXV1930
          RI     = RDIA                                                 AAXV1931
          BETA   = ALDIA                                                AAXV1932
          XK     = ZXK                                                  AAXV1933
          XVALUE(1) = 3.3D0                                             AAXV1934
          DO 20 J=2,NP                                                  AAXV1935
          XVALUE(J) = XVALUE(J-1)+.05D0                                 AAXV1936
   20     CONTINUE                                                      AAXV1937
          DO 300 II=1,NVMXP1                                            AAXV1938
          UINDEX = 2.D0*XK - 2.D0*VINDEX - 1.D0                         AAXV1939
          DO 100 IX=1,NP                                                AAXV1940
          VIBPHI(IX) = 0.D0                                             AAXV1941
          VIBLAG(IX) = 0.D0                                             AAXV1942
          XIOD   =  BETA*(XVALUE(IX)-RI)                                AAXV1943
          IF(XIOD.GT.-50.D0) THEN                                       AAXV1944
            XV   = 2.D0*XK*DEXP(-XIOD)                                  AAXV1945
            RINDEX = VINDEX                                             AAXV1946
            INDEX  = INT(RINDEX)+1                                      AAXV1947
            VIBLAG(IX) = XV**(INDEX-1)                                  AAXV1948
            COEFF  = 1.0D0                                              AAXV1949
            IF(INDEX.EQ.1)           GOTO 80                            AAXV1950
            INDXM1 = INDEX-1                                            AAXV1951
            DO  50 I=1,INDXM1                                           AAXV1952
            RA     = DFLOAT(I)                                          AAXV1953
            COEFF  =-COEFF*(VINDEX-RA+1.0D0)*(UINDEX+VINDEX-RA+1.0D0)   AAXV1954
     &                                                       /RA        AAXV1955
            VIBLAG(IX) = VIBLAG(IX)+COEFF*(XV**(VINDEX-RA))             AAXV1956
   50       CONTINUE                                                    AAXV1957
            IF(2*(INDEX/2).EQ.INDEX)     VIBLAG(IX)=-VIBLAG(IX)         AAXV1958
   80       XLGN   =-0.5D0*S14ABF(2.0D0*XK-VINDEX,IFAIL)-               AAXV1959
     -                           XV/2.0D0+(UINDEX/2.0D0)*DLOG(XV)       AAXV1960
            ALGN   = DEXP(XLGN)                                         AAXV1961
            VIBPHI(IX)=DSQRT((UINDEX*BETA)/S14AAF(VINDEX+1.D0,IFAIL))*  AAXV1962
     *                                                 ALGN*VIBLAG(IX)  AAXV1963
          END IF                                                        AAXV1964
          VIBM(IX,II) = VIBPHI(IX)                                      AAXV1965
          IF(IPRNT.GE.10)        WRITE(6,11) VIBM(IX,II),IX,II          AAXV1966
  100     CONTINUE                                                      AAXV1967
          VINDEX   = VINDEX+1.D0                                        AAXV1968
  300     CONTINUE                                                      AAXV1969
          RETURN                                                        AAXV1970
C                                                                       AAXV1971
   11     FORMAT(10X,'VIBM=',D14.6,'  IX=',I4,'  II=',I4)               AAXV1972
          END                                                           AAXV1973
C                                                                       AAXV1974
          SUBROUTINE POTD(POTRX,DPOTRX,RR,W,DW,R,UCF,WV,NBB)            AAXV1975
C*******************************************************************    AAXV1976
C USED TO GET THE DIAGONAL POTENTIAL AND DERIVATIVE                     AAXV1977
C                                                                       AAXV1978
          IMPLICIT REAL*8(A-H,O-Z)                                      AAXV1979
          PARAMETER(NDIM=96,N2DIM=4)                                    AAXV1980
          COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ        AAXV1981
          COMMON/SURF / IPNT,NPOINT,NVMXP1                              AAXV1982
          COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEPAAXV1983
     &                  NVMIN,NVMAX                                     AAXV1984
          DIMENSION RR(1),POTRX(1),DPOTRX(1)                            AAXV1985
C                                                                       AAXV1986
          CF   = UCF/(R*R)                                              AAXV1987
          NVMX = N2DIM                                                  AAXV1988
          NELM = NVMX*(NVMX+1)/2                                        AAXV1989
   90     IF(R.LT.RR(IPNT)) THEN                                        AAXV1990
            IF(IPNT.GE.2) THEN                                          AAXV1991
              IPNT = IPNT - 1                                           AAXV1992
              GOTO 90                                                   AAXV1993
            END IF                                                      AAXV1994
            WRITE(6,900) R                                              AAXV1995
            STOP 33                                                     AAXV1996
          END IF                                                        AAXV1997
  150     IF(R.GE.RR(IPNT+1)) THEN                                      AAXV1998
            IPNT = IPNT + 1                                             AAXV1999
            IF(IPNT.LT.NPOINT) GOTO 150                                 AAXV2000
            WRITE(6,910) R,RR(IPNT)                                     AAXV2001
            STOP 17                                                     AAXV2002
          END IF                                                        AAXV2003
  200     INDX   = NBB*(NBB+1)/2 + NELM*(IPNT-1)                        AAXV2004
          VMATRX = POTRX(INDX)+(R-RR(IPNT))*(POTRX(INDX+NELM)-          AAXV2005
     -                              POTRX(INDX))/(RR(IPNT+1)-RR(IPNT))  AAXV2006
          DVMTRX = DPOTRX(INDX)+(R-RR(IPNT))*(DPOTRX(INDX+NELM)-        AAXV2007
     -                             DPOTRX(INDX))/(RR(IPNT+1)-RR(IPNT))  AAXV2008
          W   = WV*WV-CF-VMATRX*BRED                                    AAXV2009
          DW  = 2.D0*CF/R-DVMTRX*BRED                                   AAXV2010
          RETURN                                                        AAXV2011
C                                                                       AAXV2012
  900     FORMAT(10X,20(1H~),' SCATTERING DISTANCE ',D14.6,'  BELOW '   AAXV2013
     &          ,'RMIN')                                                AAXV2014
  910     FORMAT(10X,20(1H~),' SCATTERING DISTANCE ',D14.7,' ABOVE '    AAXV2015
     &        ,'RMAX',D14.7)                                            AAXV2016
          END                                                           AAXV2017
C                                                                       AAXV2018
         SUBROUTINE POTND(POTRX,DPOTRX,RR,WN,DWN,R,UCF,WV,N1,N2)        AAXV2019
C*******************************************************************    AAXV2020
C USED TO GET NON-DIAGONAL MATRIX ELEMENTS AND DERIVATIVES              AAXV2021
C                                                                       AAXV2022
          IMPLICIT REAL*8(A-H,O-Z)                                      AAXV2023
          PARAMETER(NDIM=96,N2DIM=4)                                    AAXV2024
          COMMON/QUANT/ BRED,ERED,VIBQUAN,JTOT,IEDW,LWKBJ,STWKBJ        AAXV2025
          COMMON/SURF / IPNT,NPOINT,NVMXP1                              AAXV2026
          COMMON/PASS / ETOT,ESTEP,NETOT,CONV,CONVP,JJTMAX,JJTMIN,JTSTEPAAXV2027
     &                  NVMIN,NVMAX                                     AAXV2028
          DIMENSION RR(1),POTRX(1),DPOTRX(1)                            AAXV2029
C                                                                       AAXV2030
          NVMX  = N2DIM                                                 AAXV2031
          NELM  = NVMX*(NVMX+1)/2                                       AAXV2032
   90     IF(R.LT.RR(IPNT)) THEN                                        AAXV2033
            IF(IPNT.LT.2)         STOP 34                               AAXV2034
            IPNT = IPNT-1                                               AAXV2035
            GOTO 90                                                     AAXV2036
          END IF                                                        AAXV2037
  150     IF(R.GE.RR(IPNT+1)) THEN                                      AAXV2038
            IPNT = IPNT+1                                               AAXV2039
            IF(IPNT.LT.NPOINT)    GOTO 150                              AAXV2040
            STOP 18                                                     AAXV2041
          END IF                                                        AAXV2042
  200     INDX   = N1*(N1-1)/2 + N2 + NELM*(IPNT-1)                     AAXV2043
          VMATRX = POTRX(INDX)+(R-RR(IPNT))*(POTRX(INDX+NELM)-          AAXV2044
     -                    POTRX(INDX))/(RR(IPNT+1)-RR(IPNT))            AAXV2045
          DVMTRX = DPOTRX(INDX)+(R-RR(IPNT))*(DPOTRX(INDX+NELM)-        AAXV2046
     -                    DPOTRX(INDX))/(RR(IPNT+1)-RR(IPNT))           AAXV2047
          WN  = VMATRX*BRED                                             AAXV2048
          DWN = DVMTRX*BRED                                             AAXV2049
          RETURN                                                        AAXV2050
          END                                                           AAXV2051
 HE-I2   VEDW     (UNITS A.U.)                                          AAXV2052
 $EDW   RMASS=7183.033219D0, JTSTEP=30, IEDW=-1,                        AAXV2053
        EMIN=6.0D-3, ETOTAL=12.0D-3, NESTEP=0,                          AAXV2054
        IPRNT=1, NVMAX=3, VIBQUAN=9.71449D-4, MAXSTEP=100,              AAXV2055
        RVDW=7.5589075D0, DVDW=6.2473D-5, ALVDW=6.244288D-1,            AAXV2056
        RDIA=5.03927041D0, ALDIA=0.9875114D0, DDIA=5.665588D-2,         AAXV2057
        DIMASS=115691.3734D0, IOLD=0                                    AAXV2058
 $END                                                                   AAXV2059
                                                                        AAXV****
