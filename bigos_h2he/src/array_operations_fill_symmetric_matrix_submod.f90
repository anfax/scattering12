!!!*************************************************************
! 文件/File: array_operations_fill_symmetric_matrix_submod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: array_operations_fill_symmetric_matrix_submod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

submodule (array_operations_mod) array_operations_fill_symmetric_matrix_submod
   !! a submodule for append subroutines
   use utility_functions_mod, only: integer_to_character, write_message,       &
      write_error, to_lowercase
   implicit none
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine fill_symmetric_matrix_int32(matrix_, upper_lower_)
          !! Fill the upper/lower triangle of a symmetric matrix (integer).
          integer, intent(inout) :: matrix_(:,:)
          character(len = 1), intent(in) :: upper_lower_
          !--------------------------------------------------------------------!
          integer :: i_, size_1_, size_2_, size_
          !--------------------------------------------------------------------!
          size_1_ = size(matrix_, dim = 1)
          size_2_ = size(matrix_, dim = 2)
          if (size_1_ .eq. size_2_) then
             size_ = size_1_
          else
             call write_message("Error in fill_symmetric_matrix_int: size "    &
               // "in dim = 1 ("//trim(adjustl(integer_to_character(size_1_))) &
               // ") is different than in dim = 2 (" //                        &
               trim(adjustl(integer_to_character(size_2_))) // ")")
             call write_error("Adapt this subroutine to rectangle matrices")
          endif
          !--------------------------------------------------------------------!
          select case(to_lowercase(upper_lower_))
              case('l')
                  do i_ = 1, size_ - 1
                      matrix_(i_ + 1:size_, i_) = matrix_(i_, i_ + 1:size_)
                  enddo
              case('u')
                  do i_ = 1, size_ - 1
                      matrix_(i_, i_ + 1:size_) = matrix_(i_ + 1:size_, i_)
                  enddo
              case default
                  call write_message("Error: Invalid argument in " //          &
                     "fill_symmetric_matrix_int32 subroutine (upper_lower_):"  &
                     // upper_lower_)
                  call write_error("'u' or 'l' expected")
          end select
          !--------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_int32
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine fill_symmetric_matrix_sp(matrix_, upper_lower_)
         !! Fill the upper/lower triangle of a symmetric matrix
         !! (single precision).
         real(sp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
         !---------------------------------------------------------------------!
         integer :: i_, size_1_, size_2_, size_
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
             call write_message("Error in fill_symmetric_matrix_sp: size "     &
               // "in dim = 1 ("//trim(adjustl(integer_to_character(size_1_))) &
               // ") is different than in dim = 2 (" //                        &
               trim(adjustl(integer_to_character(size_2_))) // ")")
             call write_error("Adapt this subroutine to rectangle matrices")
         endif
         !---------------------------------------------------------------------!
         select case(to_lowercase(upper_lower_))
             case('l')
                 do i_ = 1, size_ - 1
                     matrix_(i_ + 1:size_, i_) = matrix_(i_, i_ + 1:size_)
                 enddo
             case('u')
                 do i_ = 1, size_ - 1
                     matrix_(i_, i_ + 1:size_) = matrix_(i_ + 1:size_, i_)
                 enddo
             case default
                  call write_message("Error: Invalid argument in " //          &
                     "fill_symmetric_matrix_int32 subroutine (upper_lower_):"  &
                     // upper_lower_)
                  call write_error("'u' or 'l' expected")
         end select
         !---------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_sp
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      module subroutine fill_symmetric_matrix_dp(matrix_, upper_lower_)
         !! fill the upper/lower triangle of a symmetric matrix
         real(dp), intent(inout) :: matrix_(:,:)
         character(len = 1), intent(in) :: upper_lower_
         !---------------------------------------------------------------------!
         integer :: i_, size_1_, size_2_, size_
         !---------------------------------------------------------------------!
         size_1_ = size(matrix_, dim = 1)
         size_2_ = size(matrix_, dim = 2)
         if (size_1_ .eq. size_2_) then
            size_ = size_1_
         else
             call write_message("Error in fill_symmetric_matrix_dp: size "     &
               // "in dim = 1 ("//trim(adjustl(integer_to_character(size_1_))) &
               // ") is different than in dim = 2 (" //                        &
               trim(adjustl(integer_to_character(size_2_))) // ")")
             call write_error("Adapt this subroutine to rectangle matrices")
         endif
         !---------------------------------------------------------------------!
         select case(to_lowercase(upper_lower_))
            case('l')
             do i_ = 1, size_ - 1
                matrix_(i_ + 1 : size_, i_) = matrix_(i_, i_ + 1 : size_)
             enddo
            case('u')
             do i_ = 1, size_ - 1
                matrix_(i_, i_ + 1: size_) = matrix_(i_ + 1 : size_, i_)
             enddo
            case default
                  call write_message("Error: Invalid argument in " //          &
                     "fill_symmetric_matrix_int32 subroutine (upper_lower_):"  &
                     // upper_lower_)
                  call write_error("'u' or 'l' expected")
         end select
         !---------------------------------------------------------------------!
      end subroutine fill_symmetric_matrix_dp
   !---------------------------------------------------------------------------!
end submodule array_operations_fill_symmetric_matrix_submod
