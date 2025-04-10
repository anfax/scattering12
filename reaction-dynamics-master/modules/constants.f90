!!!*************************************************************
! 文件/File: constants.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: constants.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

MODULE constants
    ! States and gridpoints
    integer*8,parameter::vstate=0
    integer*8,parameter::lstate=0
    integer*8,parameter::mstate=0
    integer*8,parameter::nR = 100
    integer*8,parameter::nr1 = 100
    integer*8,parameter::ntheta = 20
    integer*8,parameter::ncheb = 1

    real*8,parameter::r1min=0.5
    real*8,parameter::r1max=10.5
    real*8,parameter::Rmin=0.5
    real*8,parameter::Rmax=10.5
    real*8,parameter::delR=(10.5-0.5)/(nR-1.0)
    real*8,parameter::delr1=(10.5-0.5)/(nr1-1.0)
    real*8,parameter::Ne=1.0
    real*8,parameter::Ro=6.0
    real*8,parameter::dell=0.15
    real*8,parameter::ko=35.85
    real*8,parameter::hcross=1.0
    real*8,parameter::Hmax = 1.5
    real*8,parameter::Hmin = 0.0

    ! All the masses and stuff
    real*8,parameter::amass=1.007
    real*8,parameter::bmass=1.007
    real*8,parameter::cmass=1.007
    real*8,parameter::amu = 1822.89d0
    real*8,parameter::rmR = amu*amass*(bmass+cmass)/(amass+bmass+cmass)
    real*8,parameter::rmr1 = amu*bmass*cmass/(bmass+cmass)
    real*8,parameter::rm = amu*bmass*cmass/(bmass+cmass)
END MODULE constants
