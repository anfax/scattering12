/**************************************************************
 * 文件/File: dipolsys.h
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: dipolsys.h
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:42
 **************************************************************
 */

c ---- Quantum number parameters for stretching wave functions ----- 

c -----  Morse parameter Kappa  = sqrt(2*m*D)/a/hbar (see Spirko, JMS, 112, 183(1985))
      REAL*8  KAPPA(2,3),RMin(2),RMax(2),RhoMaximum,
     1        StrengthLimit,AbsTemper,TempCO,IntensLimit,
     2        EnerMaxAB,EnerMaxC,NuMin,Numax
      DOUBLE PRECISION Te_A,Te_C,ZeroEnergyA,ZeroEnergyC

      INTEGER NV1V3Max(3),Nvib(3,2),NumPtsBond,
     2        MaxNumberV1V3,MaxV1V3,
     1        MinNVibAB,MaxNVibAB,
     2        MinNVibC,MaxNVibC, 
     3        MinGamVibAB,MaxGamVibAB,
     4        MinGamVibC,MaxGamVibC,
     5        MinV2AB,MaxV2AB,
     6        MinV2C,MaxV2C,MinSurf,MaxSurf,
     7        MinJAB,MaxJAB,
     8        MinJC,MaxJC,
     7        MinKAB,MaxKAB,MinKC,
     8        MaxKC,CalcOrRead,MaxNvibMax,OUT_COEFF,
     9        NumFileLetters


      INTEGER KMAXAB,KMAXC,NV2MaxAB,NV2MaxC,
     1        NumPtsRho,NumDirLetters,NV2Max(2)

      COMMON /IntDIPOL/ NV1V3Max,Nvib,
     1        NumPtsBond,MaxNumberV1V3,
     2        MaxV1V3,KMAXAB,KMAXC,NV2MaxAB,NV2MaxC,
     3        NumPtsRho,
     1        MinNVibAB,MaxNVibAB,
     2        MinNVibC,MaxNVibC, 
     3        MinGamVibAB,MaxGamVibAB,
     4        MinGamVibC,MaxGamVibC,
     5        MinV2AB,MaxV2AB,
     6        MinV2C,MaxV2C,MinSurf,MaxSurf,
     7        MinJAB,MaxJAB,
     8        MinJC,MaxJC,
     7        MinKAB,MaxKAB,MinKC,
     8        MaxKC,CalcOrRead,MaxNvibMax,NV2Max,OUT_COEFF,
     9        NumFileLetters

      COMMON /REALDIP/ KAPPA,RMin,RMax,RhoMaximum,
     1        AbsTemper,TempCO,IntensLimit,
     2        StrengthLimit,Te_A,Te_C,ZeroEnergyA,ZeroEnergyC,
     2        EnerMaxAB,EnerMaxC,NuMin,Numax,NumDirLetters


      CHARACTER*80 DirInput,FileOutput
      COMMON /STRING/ DirInput,FileOutput

    
