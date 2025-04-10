!!!*************************************************************
! 文件/File: syusr.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: syusr.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck syusr
#include "hiutil.inc.F90"
subroutine syusr (irpot, readpt, iread)
!  dummy syusr subroutine
implicit none
integer, intent(inout) :: irpot
logical, intent(inout) :: readpt
integer, intent(in) :: iread
character*(*), intent(in) :: fname
UNUSED(irpot)
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
entry ptrusr (fname, readpt)
UNUSED(fname)
entry savusr (readpt)
entry chkusr
return
end
