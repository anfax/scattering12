!!!*************************************************************
! 文件/File: grovec.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: grovec.F90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

#ifdef HAS_CONCATENATION_OPERATOR
#define CONCATENATE(a,b) a ## b
#else
#define IDENTITY(x) x
#define CONCATENATE(a,b) IDENTITY(a)IDENTITY(b)
#endif

! grovec : GROwable VECtor
! grovec provides a way to store a sparse array of size n in the form of blocks that are allocated when needed
module mod_grovec
  implicit none
  private
  integer, public, protected :: g_num_allocated_blocks = 0
  public :: print_grovec_stats

#define GROVEC_CLASS_NAME dgrovec_type
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_type.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

#define GROVEC_CLASS_NAME igrovec_type
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_type.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

contains

#define GROVEC_CLASS_NAME dgrovec_type
#define GROVEC_ELEMENT_TYPE real(8)
#include "grovec_imp.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

#define GROVEC_CLASS_NAME igrovec_type
#define GROVEC_ELEMENT_TYPE integer
#include "grovec_imp.F90"
#undef GROVEC_ELEMENT_TYPE
#undef GROVEC_CLASS_NAME

subroutine print_grovec_stats(unit)
  integer, intent(in) :: unit
  write(unit, *) 'number of allocated grovec blocks: ', g_num_allocated_blocks
end subroutine



end module mod_grovec
