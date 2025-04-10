!!!*************************************************************
! 文件/File: struct_mat.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: struct_mat.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine struct_mat
  use Brep, only : p_Patch 
  use Matrices, only : gamma_c,gamma
  implicit none 
  integer,parameter :: dime=216
  integer i,j
  do i=1,p_Patch
     do j=1,dime
        if (gamma(i,j).ge.1.0d-8) then
           !write(318,*)i,gamma_c(i,j)
        end if
     end do
  end do

end subroutine struct_mat
