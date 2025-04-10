!!!*************************************************************
! 文件/File: main.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: main.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


program main 
  !use bsp_info
  implicit none
  !call bspbegin(bspnprocs())

  ! Main program of code FERM3D.x 
write(6,*)'start FERM3D code'
  call grid ! grid calculation


  call over_calc ! surface harmonics calculation

  call potential ! potential setup

  call mat_el ! matrix elements and solution of R-matrix

write(6,*)'end FERM3D code'
  !call bspend()


end program main
