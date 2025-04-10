!!!*************************************************************
! 文件/File: efield.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: efield.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


REAL(kind=8) FUNCTION efield(tau)
  implicit none

  INTEGER :: cycles
  REAL(kind=8) :: tau, omega, t, pi

  pi = 4.d0*atan(1.d0)
  omega = 0.55d0
  cycles = 50
  t = cycles*pi/omega

  efield = ((sin(0.5d0*pi*tau/T))**2)*sin(omega*tau)
  RETURN
END FUNCTION efield
