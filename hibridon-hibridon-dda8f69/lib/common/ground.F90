!!!*************************************************************
! 文件/File: ground.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: ground.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

!comdeck ground
#include "hiutil.inc.F90"
subroutine ground(wf, r, nch, nphoto, mxphot)
implicit none
real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
! ground state wavefunction in each of nch components
! nphoto is number of difference ground state wavefunctions
real(8), intent(out) :: r  ! value of separation coordinate
integer, intent(in) :: nch  ! total number of channels (row dimension of q)
integer, intent(in) :: nphoto ! number of different wavefunctions calculated
! column index of q vector
integer, intent(in) :: mxphot  ! maximum size of q vector (mxphot .ge. nch*nphoto)

UNUSED(wf)
UNUSED(r)
UNUSED(nch)
UNUSED(nphoto)
UNUSED(mxphot)
end subroutine

subroutine wfintern(wf, yymin, nch, nphoto, nny, ifull)
implicit none
real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
! ground state wavefunction in each of nch components
! nphoto is number of difference ground state wavefunctions
real(8), intent(in) :: yymin
integer, intent(in) :: nch  ! total number of channels (row dimension of q)
integer, intent(in) :: nphoto ! number of different wavefunctions calculated
! column index of q vector
integer, intent(in) :: nny
logical, intent(in) :: ifull

UNUSED(wf)
UNUSED(yymin)
UNUSED(nch)
UNUSED(nny)
UNUSED(ifull)

return
end subroutine
