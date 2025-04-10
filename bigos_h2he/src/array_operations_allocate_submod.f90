!!!*************************************************************
! 文件/File: array_operations_allocate_submod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: array_operations_allocate_submod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

submodule (array_operations_mod) array_operations_allocate_submod
   !! a submodule for allocate subroutines
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                              1-D case
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_1d_int32(array_, size_)
         !! allocate a 1d array and fill it with 0s (intger version)
         integer(int32), allocatable, intent(inout) :: array_(:)
         integer(int32), intent(in)                 :: size_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_1d_int32
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_1d_sp(array_, size_)
         !! allocate a 1d array and fill it with 0s (single precision version)
         real(sp), allocatable, intent(inout)  :: array_(:)
         integer(int32), intent(in)            :: size_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_1d_sp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_1d_dp(array_, size_)
         !! allocate a 1d array and fill it with 0s (double precision version)
         real(dp), allocatable, intent(inout)  :: array_(:)
         integer(int32), intent(in)            :: size_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_1d_dp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                              2-D case
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_2d_int32(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (intger version)
         integer(int32), allocatable, intent(inout) :: array_(:,:)
         integer(int32), intent(in)                 :: size1_, size2_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_2d_int32
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_2d_sp(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (single precision version)
         real(sp), allocatable, intent(inout)  :: array_(:,:)
         integer(int32), intent(in)            :: size1_, size2_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_2d_sp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_2d_dp(array_, size1_, size2_)
         !! allocate a 2d array and fill it with 0s (double precision version)
         real(dp), allocatable, intent(inout)  :: array_(:,:)
         integer(int32), intent(in)            :: size1_, size2_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_2d_dp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                              3-D case
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_3d_int32(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (intger version)
         integer(int32), allocatable, intent(inout) :: array_(:,:,:)
         integer(int32), intent(in)                 :: size1_, size2_, size3_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_, size3_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_3d_int32
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_3d_sp(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (single precision version)
         real(sp), allocatable, intent(inout)  :: array_(:,:,:)
         integer(int32), intent(in)            :: size1_, size2_, size3_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_, size3_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_3d_sp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      pure module subroutine allocate_3d_dp(array_, size1_, size2_, size3_)
         !! allocate a 3d array and fill it with 0s (double precision version)
         real(dp), allocatable, intent(inout)  :: array_(:,:,:)
         integer(int32), intent(in)            :: size1_, size2_, size3_
         !---------------------------------------------------------------------!
         if (allocated(array_)) deallocate(array_)
         allocate(array_(size1_, size2_, size3_))
         array_ = 0
         !---------------------------------------------------------------------!
      end subroutine allocate_3d_dp
      !------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
end submodule array_operations_allocate_submod
