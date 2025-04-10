!!!*************************************************************
! 文件/File: assert.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: assert.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine fassert(filename,linenum)
	character(*), intent(in) :: filename
	integer,      intent(in) :: linenum
	write(6,'(a,a,a,i4.4)') "assertion failed in file ", filename,":",linenum
	stop 1
end subroutine fassert

! subroutine FortranAssert(prepost,expression,filename,linenum)
! 	character(*), intent(in) :: prepost,expression, filename
! 	integer,      intent(in) :: linenum
! 	write(6,'(a,a,a,a,a,a,i4.4)') prepost," (", expression, ") failed in file ", filename," line ",linenum
! 	stop
! end subroutine FortranAssert
