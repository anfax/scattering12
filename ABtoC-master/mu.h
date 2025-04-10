/**************************************************************
 * 文件/File: mu.h
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: mu.h
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:42
 **************************************************************
 */

      INTEGER MaxNumVibMatEl,MaxNumCoeffAB,MaxNumCoeffC,
     1        MaxNumRhoPoints,MaxNumLevelAB,MaxNumLevelC,
     2        UpperV2,UpperK,UpperNvib

      INTEGER NRECUN
      PARAMETER ( NRECUN = 8 )
      PARAMETER(MaxNumVibMatEl=9100000)
      PARAMETER(MaxNumCoeffAB=40000)
      PARAMETER(MaxNumCoeffC =20000)
      PARAMETER(MaxNumLevelAB=30000)
      PARAMETER(MaxNumLevelC =40000)

      PARAMETER(UpperV2   = 31)
      PARAMETER(UpperK    = 11)
      PARAMETER(UpperNvib = 35)

      REAL*8 AlmostZero
      PARAMETER(AlmostZero=1D-16)

      INTEGER InpUnit,OutUnit,TMPUnit,MaxNumJ
      
      
      PARAMETER(MaxNumRhoPoints=2501)
      PARAMETER(MaxNumJ=20)


      INTEGER NotateStretch(3,2,UpperNvib)  

      INTEGER :: Gns(4) ! = (/1,1,3,3/)

      COMMON /IntMu/ InpUnit,OutUnit,TMPUnit,
     1        NotateStretch,gns,GammaLevel

      CHARACTER*2 GammaLevel(4)

        LOGICAL FlagCheck 
        PARAMETER(FlagCheck = .FALSE.)

      CHARACTER*20 FileStretchAB,
     1            FileStretchC,FileGamma,FileBendAB,
     2            FileBendC,FileCoefAB,FileCoefC,FileTMP,
     2            FileVibME

      DATA 
     2    FileStretchAB /'stret_ab.dat'/,
     3    FileStretchC  /'stret_c.dat'/,
     4    FileGamma     /'gamma_ab.dat'/,
     5    FileBendAB    /'bend_ab.dat'/,
     6    FileBendC     /'bend_c.dat'/,
     7    FileCoefAB    /'coef_ab.dat'/,
     8    FileCoefC     /'coef_c.dat'/,
     9    FileVibME     /'vib_me.tmp'/,
     1    FileTMP       /'intens.dat'/


       INTEGER NFileMu  
      PARAMETER( NFileMu=80 )

      INTEGER StartTime,TimePassed
      COMMON /CURRTIME/ StartTime,TimePassed




