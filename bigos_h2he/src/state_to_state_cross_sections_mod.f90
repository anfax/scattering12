!!!*************************************************************
! 文件/File: state_to_state_cross_sections_mod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: state_to_state_cross_sections_mod.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

module state_to_state_cross_sections_mod
   !! This module provides functions and subroutines for calculating and
   !! analyzing state-to-state cross-sections. It is divided into three
   !! parts:
   !---------------------------------------------------------------------------!
   !! (1) Calculating cross-sections: Functions for computing state-to-state
   !!    cross-sections based on quantum states, S-matrix, and scattering
   !!    parameters ("calculate_state_to_state_cross_section",
   !!    "compute_individual_cross_section", "get_block_indices",
   !!    "sum_cross_section_contributions", "compute_real_component",
   !!    "compute_imag_component")
   !! (2) Printing cross-sections: Subroutines to output the largest partial
   !!    cross-sections, providing both basic and detailed information.
   !!    ("print_largest_partial_cross_sections",
   !!    "print_basic_cross_section_info", "print_detailed_cross_section_info")
   !! (3) Threshold checking: Subroutine to check if the computed cross-sections
   !!    meet specified convergence conditions
   !!    ("check_cross_section_thresholds")
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   use utility_functions_mod, only: write_message, time_count_summary,         &
      integer_to_character, float_to_character, file_io_status
   use array_operations_mod, only: allocate_1d
   use global_variables_mod
   use physics_utilities_mod, only: total_energy
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   private
   public :: initialize_cross_section_arrays,                                  &
      calculate_state_to_state_cross_section, add_cross_sections,              &
      print_largest_partial_cross_sections, print_cross_sections_for_jtot,     &
      print_final_cross_sections, save_partial_xs_file_header,                 &
      save_partial_xs_single_block, determine_largest_cross_sections,          &
      check_cross_section_thresholds
   !---------------------------------------------------------------------------!
   contains
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                   Initialization of cross-sections arrays
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine initialize_cross_section_arrays(size_, partial_block_,        &
         partial_jtot_, accumulated_)
         !! allocate arrays keeping accumulated and partial cross-sections
         !! in each jtot and parity block
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: size_
            !! number of open basis levels (defines the size of both arrays)
         real(dp), intent(inout), allocatable :: partial_block_(:)
            !! array holding partial cross-sections in a parity block
         real(dp), intent(inout), allocatable :: partial_jtot_(:)
            !! array holding partial cross-sections in a jtot block
         real(dp), intent(inout), allocatable :: accumulated_(:)
            !! array holding accumulated cross-sections
         !---------------------------------------------------------------------!
         call allocate_1d(partial_block_, size_*size_)
         call allocate_1d(partial_jtot_, size_*size_)
         call allocate_1d(accumulated_, size_*size_)
         !---------------------------------------------------------------------!
      end subroutine initialize_cross_section_arrays
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                         Calculating cross-sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine calculate_state_to_state_cross_section(                       &
         total_angular_momentum_, open_basis_levels_, basis_wavevectors_,      &
         s_matrix_real_, s_matrix_imag_, channel_indices_, channel_l_values_,  &
         cross_section_array_)
         !! Calculates all state-to-state cross-sections.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! holds indices to the basis arrays which correspond to open channels 
         real(dp), intent(in) :: basis_wavevectors_(:)
            !! holds wavevectors k_{i}
         real(dp), intent(in) :: s_matrix_real_(:,:), s_matrix_imag_(:,:)
            !! real and imaginary parts of the S-matrix
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channel_l_values_(:)
            !! holds all values of l
         real(dp), intent(inout) :: cross_section_array_(:)
            !! array holding all XSs
         !---------------------------------------------------------------------!
         integer(int32) :: number_of_open_basis_levels_, initial_state_,       &
            final_state_, cross_section_index_
         real(dp) :: start_time, finish_time, calculation_time
         !---------------------------------------------------------------------!
         call CPU_TIME(start_time)
         number_of_open_basis_levels_ = size(open_basis_levels_)
         cross_section_array_ = 0
         !---------------------------------------------------------------------!
         do initial_state_ = 1, number_of_open_basis_levels_
            do final_state_ = 1, number_of_open_basis_levels_                   
               cross_section_index_ = (initial_state_-1)                       &
                  * number_of_open_basis_levels_ + final_state_
               cross_section_array_(cross_section_index_) =                    &
                  compute_individual_cross_section(initial_state_,final_state_,&
                     open_basis_levels_, basis_wavevectors_, s_matrix_real_,   &
                     s_matrix_imag_, channel_indices_, channel_l_values_,      &
                     total_angular_momentum_)
            enddo
         enddo
         !---------------------------------------------------------------------!
         CALL CPU_TIME(finish_time)
         if (print_level >= 2) then
            call time_count_summary(start_time, finish_time, calculation_time, &
               "Cross-sections calculations completed in ")
         endif
         !---------------------------------------------------------------------!
      end
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function compute_individual_cross_section(initial_state_, final_state_,  &
         open_basis_levels_, basis_wavevectors_, s_matrix_real_,          &
         s_matrix_imag_, channel_indices_, channel_l_values_,                  &
         total_angular_momentum_) result(cross_section_)
         !! Calculates cross-section for a given initial and final state.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: initial_state_
            !! index pointing to the intial state in basis arrays
         integer(int32), intent(in) :: final_state_
            !! index pointing to the final state in basis arrays
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! holds indices to the basis arrays which correspond to open channels
         real(dp), intent(in) :: basis_wavevectors_(:)
            !! holds wavevectors k_{i}
         real(dp), intent(in) :: s_matrix_real_(:,:), s_matrix_imag_(:,:)
            !! real and imaginary parts of the S-matrix
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), intent(in) :: channel_l_values_(:)
            !! holds all values of l
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp) :: cross_section_
            !! (output) cross-section
         !---------------------------------------------------------------------!
         integer(int32) :: initial_index_, final_index_, v_initial_,           &
            j_initial_, v_final_, j_final_
         real(dp) :: wavevector_initial_, sum_contributions_
         integer(int32), allocatable :: init_block_indices_(:),                &
            final_block_indices_(:)
         !---------------------------------------------------------------------!
         initial_index_ = open_basis_levels_(initial_state_)
         v_initial_ = vib_levels(initial_index_)
         j_initial_ = rot_levels(initial_index_)
         wavevector_initial_ = basis_wavevectors_(initial_state_)
         !---------------------------------------------------------------------!
         final_index_ = open_basis_levels_(final_state_)
         v_final_ = vib_levels(final_index_)
         j_final_ = rot_levels(final_index_)
         !---------------------------------------------------------------------!
         init_block_indices_ = get_block_indices(v_initial_, j_initial_,       &
            channel_indices_)
         final_block_indices_ = get_block_indices(v_final_, j_final_,          &
            channel_indices_)
         !---------------------------------------------------------------------!
         sum_contributions_ = sum_cross_section_contributions(                 &
            init_block_indices_, final_block_indices_, s_matrix_real_,         &
            s_matrix_imag_, channel_l_values_)
         !---------------------------------------------------------------------!
         cross_section_ = sum_contributions_                                   &
            * real(2 * total_angular_momentum_ + 1, dp)                        &
            / (real(2 * j_initial_ + 1, dp) * wavevector_initial_**2.0_dp) * pi
         !---------------------------------------------------------------------!
      end function compute_individual_cross_section
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function get_block_indices(v_, j_, channel_indices_) result(block_indices_)
         !! Returns indices from channel_indices_ that match the specified quantum state.
         integer(int32), intent(in) :: v_
            !! vibrational quantum number
         integer(int32), intent(in) :: j_
            !! rotational quantum number
         integer(int32), intent(in) :: channel_indices_(:)
            !! holds the indices pointing to the basis arrays
         integer(int32), allocatable :: block_indices_(:)
            !! (output) indices that match the specified quantum state
         !---------------------------------------------------------------------!
         integer(int32) :: number_of_channels_, match_count_, index_
         !---------------------------------------------------------------------!
         number_of_channels_ = size(channel_indices_)
         match_count_ = 0
         !---------------------------------------------------------------------!
         ! Count the number of matches to preallocate the array
         !---------------------------------------------------------------------!
         do index_ = 1, number_of_channels_
            if (vib_levels(channel_indices_(index_)) == v_ .and.               &
               rot_levels(channel_indices_(index_)) == j_) then
               match_count_ = match_count_ + 1
            endif
         enddo
         !---------------------------------------------------------------------!
         ! Allocate array with the correct size
         !---------------------------------------------------------------------!
         allocate(block_indices_(match_count_))
         match_count_ = 0
         !---------------------------------------------------------------------!
         ! Fill the array with matching indices.
         !---------------------------------------------------------------------!
         do index_ = 1, number_of_channels_
            if (vib_levels(channel_indices_(index_)) == v_ .and.               &
               rot_levels(channel_indices_(index_)) == j_) then
               match_count_ = match_count_ + 1
               block_indices_(match_count_) = index_
            endif
         enddo
         !---------------------------------------------------------------------!
      end function get_block_indices
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function sum_cross_section_contributions(init_indices_, final_indices_,  &
         s_matrix_real_, s_matrix_imag_, channel_l_values_) result(sum_contributions_)
         !! Sum the contributions to the cross-section from S-matrix components
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: init_indices_(:)
            !! indices pointing to basis element involving initial state
         integer(int32), intent(in) :: final_indices_(:)
            !! indices pointing to basis element involving final state
         integer(int32), intent(in) :: channel_l_values_(:)
            !! holds all values of l
         real(dp), intent(in) :: s_matrix_real_(:,:), s_matrix_imag_(:,:)
            !! real and imaginary parts of the S-matrix
         real(dp) :: sum_contributions_
            !! (output) contribution to the cross-section from S-matrix components
         !---------------------------------------------------------------------!
         integer(int32) :: initial_index_, final_index_, l_initial_, l_final_
         real(dp) :: term_real_, term_imag_, term_squared_
         !---------------------------------------------------------------------!
         sum_contributions_ = 0.0_dp

         do initial_index_ = 1, size(init_indices_)
           l_initial_ = channel_l_values_(init_indices_(initial_index_))
           do final_index_ = 1, size(final_indices_)
               l_final_ = channel_l_values_(final_indices_(final_index_))

               term_real_ = compute_real_component(                            &
                  init_indices_(initial_index_), final_indices_(final_index_), &
                  l_initial_, l_final_, s_matrix_real_)
               term_imag_ = compute_imag_component(                            &
                  init_indices_(initial_index_), final_indices_(final_index_), &
                  s_matrix_imag_)

               term_squared_ = term_real_**2.0_dp + term_imag_**2.0_dp
               sum_contributions_ = sum_contributions_ + term_squared_
           end do
         end do
      end function sum_cross_section_contributions
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function compute_real_component(initial_index_, final_index_, l_initial_,&
         l_final_, s_matrix_real_) result(real_contribution_)
         !! Computes the real part of the cross-section contribution
         !! for given indices.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: initial_index_
            !! index pointing to the basis element involving initial state
         integer(int32), intent(in) :: final_index_
            !! index pointing to the basis element involving final state
         integer(int32), intent(in) :: l_initial_
            !! pre-collisional \\(l\\) value
         integer(int32), intent(in) :: l_final_
            !! post-collisional \\(l\\) value
         real(dp), intent(in) :: s_matrix_real_(:,:)
            !! real part of the S-matrix
         real(dp) :: real_contribution_
            !! (output) contribution to the cross-section from real part
            !! of the S-matrix
         !---------------------------------------------------------------------!
         if (l_initial_ == l_final_) then
           if (initial_index_ == final_index_) then
               real_contribution_ = 1.0_dp - s_matrix_real_(initial_index_, final_index_)
           else
               real_contribution_ = - s_matrix_real_(initial_index_, final_index_)
           endif
         else
           real_contribution_ = - s_matrix_real_(initial_index_, final_index_)
         endif
      end function compute_real_component
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      function compute_imag_component(initial_index_, final_index_,            &
         s_matrix_imag_) result(imag_contribution_)
         !! Computes the imaginary part of the cross-section contribution
         !! for given indices.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: initial_index_
            !! index pointing to the basis element involving initial state
         integer(int32), intent(in) :: final_index_
            !! index pointing to the basis element involving final state
         real(dp), intent(in) :: s_matrix_imag_(:,:)
            !! imaginary part of the S-matrix
         real(dp) :: imag_contribution_
            !! (output) contribution to the cross-section from imaginary
            !! part of the S-matrix
         !---------------------------------------------------------------------!
         imag_contribution_ = -s_matrix_imag_(initial_index_, final_index_)
         !---------------------------------------------------------------------!
      end function compute_imag_component
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                         Adding cross-sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine add_cross_sections(size_, partial_, accumulated_)
         !! Add partial cross-sections to accumulated cross-sections
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: size_
            !! number of open basis levels (defines the size of both arrays)
         real(dp), intent(in) :: partial_(size_*size_)
            !! array holding partial cross-sections
         real(dp), intent(inout) :: accumulated_(size_*size_)
            !! array holding accumulated cross-sections
         !---------------------------------------------------------------------!
         integer(int32) :: index_1_, index_2_, cross_section_index_
         !---------------------------------------------------------------------!
         do index_1_ = 1, size_
            do index_2_ = 1, size_
               cross_section_index_ = (index_1_-1) * size_ + index_2_
               accumulated_(cross_section_index_) =                            &
                  accumulated_(cross_section_index_)                           &
                  + partial_(cross_section_index_)
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine add_cross_sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!    
      !                        Printing cross-sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_largest_partial_cross_sections(total_angular_momentum_, &
         largest_elastic_xs_, largest_inelastic_xs_, elastic_index_,           &
         inelastic_index_1_, inelastic_index_2_, open_basis_levels_)
         !! Print the largest partial elastic and inelastic state-to-state
         !! cross-sections in a given block.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(in) :: largest_elastic_xs_
            !! the largest partial elastic state-to-state XS in the block
         real(dp), intent(in) :: largest_inelastic_xs_
            !! the largest partial inelastic state-to-state XS in the block
         integer(int32), intent(in) :: elastic_index_
            !! index pointing indirectly to quantum numbers associated with
            !! the largest partial elastic state-to-state XS in the block
         integer(int32), intent(in) :: inelastic_index_1_, inelastic_index_2_
            !! indices pointing indirectly to quantum numbers associated with
            !! the largest partial inelastic state-to-state XS in the block
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels 
         !---------------------------------------------------------------------!
         if (print_level <= 2) then
            call print_basic_cross_section_info(total_angular_momentum_,       &
               largest_elastic_xs_, "elastic")
            call print_basic_cross_section_info(total_angular_momentum_,       &
               largest_inelastic_xs_, "inelastic")
         else if (print_level >= 3) then
            call print_detailed_cross_section_info(total_angular_momentum_,    &
               largest_elastic_xs_, elastic_index_, elastic_index_,            &
               open_basis_levels_, "elastic")
            call print_detailed_cross_section_info(total_angular_momentum_,    &
               largest_inelastic_xs_, inelastic_index_1_, inelastic_index_2_,  &
               open_basis_levels_, "inelastic")
         endif
         !---------------------------------------------------------------------!
      end subroutine print_largest_partial_cross_sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_basic_cross_section_info(total_angular_momentum_,       &
         cross_section_value_, type_label_)
         !! Prints basic information about the largest elastic and inelastic
         !! state-to-state xs (print_level <= 2)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         real(dp), intent(in) :: cross_section_value_
            !! value of the cross-section
         character(len=*), intent(in) :: type_label_
            !! "elastic" or "inelastic"
         !---------------------------------------------------------------------!
         call write_message(" - Largest partial " // trim(type_label_) //      &
            " state-to-state for JTOT = " //                                   &
            trim(adjustl(integer_to_character(total_angular_momentum_))) //    &
            ": " // trim(adjustl(float_to_character(cross_section_value_))))
         !---------------------------------------------------------------------!
      end subroutine print_basic_cross_section_info
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_detailed_cross_section_info(total_angular_momentum_,    &
         cross_section_value_, index_1_, index_2_, open_basis_levels_,         &
         type_label_)
         !! Prints detailed information about the largest elastic and inelastic
         !! state-to-state xs (print_level >= 3)
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: index_1_, index_2_
            !! indices pointing indirectly to quantum numbers associated with
            !! the largest partial XS in the block
         real(dp), intent(in) :: cross_section_value_
            !! value of the cross-section
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels 
         character(len=*), intent(in) :: type_label_
            !! "elastic" or "inelastic"
         !---------------------------------------------------------------------!
         character(len=200) :: header_line_, line_
         !---------------------------------------------------------------------!
         call write_message("Largest partial " // trim(type_label_) //         &
            " state-to-state for JTOT = " //                                   &
            trim(adjustl(integer_to_character(total_angular_momentum_))))
         write(header_line_, "(2x,a4,2x,a4,2x,a2,2x,a4,2x,a4,16x,a2)")         &
            "v1_f", "j1_f", "<-", "v1_i", "j1_i", "XS"
         call write_message(header_line_)
         !---------------------------------------------------------------------!
         write(line_, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8)")                    &
           vib_levels(open_basis_levels_(index_1_)),                           &
           rot_levels(open_basis_levels_(index_1_)),                           &
           vib_levels(open_basis_levels_(index_2_)),                           &
           rot_levels(open_basis_levels_(index_2_)), cross_section_value_
         call write_message(line_)
         !---------------------------------------------------------------------!
      end subroutine print_detailed_cross_section_info
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_cross_sections_for_jtot(total_angular_momentum_,        &
         open_basis_levels_, cross_sections_)
         !! Prints information about cross-sections at the end of each
         !! total angular momentum (jtot) block
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: total_angular_momentum_
            !! total angular momentum
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels
         real(dp), intent(in) :: cross_sections_(:)
            !! holds values of the cross-sections
         !---------------------------------------------------------------------!
         call write_message("Cross sections for J: "//                         &
            trim(adjustl(integer_to_character(total_angular_momentum_))) //    &
            " and energy: " //                                                 &
            trim(adjustl(float_to_character(total_energy()*hartree_to_cm,      &
            "(F10.4)"))) // " cm-1")
         !---------------------------------------------------------------------!
         call print_all_cross_sections(open_basis_levels_, cross_sections_)
         !---------------------------------------------------------------------!
      end subroutine print_cross_sections_for_jtot
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_final_cross_sections(open_basis_levels_, cross_sections_)
         !! Prints information about cross-sections at the end of the program
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels
         real(dp), intent(in) :: cross_sections_(:)
            !! holds values of the cross-sections
         !---------------------------------------------------------------------!
         call write_message("Final state-to-state XS")
         !---------------------------------------------------------------------!
         call print_all_cross_sections(open_basis_levels_, cross_sections_)
         !---------------------------------------------------------------------!
      end subroutine print_final_cross_sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine print_all_cross_sections(open_basis_levels_, cross_sections_)
         !! Prints information about cross-sections from provided
         !! "cross_sections_" array
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels
         real(dp), intent(in) :: cross_sections_(:)
            !! holds values of the cross-sections
         !---------------------------------------------------------------------!
         character(len=200) :: header_line_, line_
         integer(int32) :: index_1_, index_2_, cross_section_index_,           &
            number_of_open_basis_levels_
         !---------------------------------------------------------------------!
         write(header_line_, "(2x,a4,2x,a4,2x,a2,2x,a4,2x,a4,14x,a4,16x,a2)")  &
            "v1_f", "j1_f", "<-", "v1_i", "j1_i", "Ekin", "XS"
         call write_message(header_line_)
         !---------------------------------------------------------------------!
         number_of_open_basis_levels_ = size(open_basis_levels_)
         !---------------------------------------------------------------------!
         do index_1_ = 1, size(open_basis_levels_)
            do index_2_ = 1, size(open_basis_levels_)
               cross_section_index_ = (index_1_ - 1)                           &
                  * number_of_open_basis_levels_ + index_2_
               write(line_, "(2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8,2X,E16.8)")     &
                  vib_levels(open_basis_levels_(index_2_)),                    &
                  rot_levels(open_basis_levels_(index_2_)),                    &
                  vib_levels(open_basis_levels_(index_1_)),                    &
                  rot_levels(open_basis_levels_(index_1_)),                    &
                  (total_energy()-internal_energies(open_basis_levels_(index_1_)))&
                     * hartree_to_cm, cross_sections_(cross_section_index_)
               call write_message(line_)
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine print_all_cross_sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                        Saving partial cross-sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine save_partial_xs_file_header
         !! save "header" of the partial cross-sections file:
         !! -- label, "itype", number of levels in the basis, reduced mass of the system
         !! -- vibrational and rotational quantum numbers
         !! -- rovibrational energies
         !! -- index pointing to the initial level and the kinetic/total energy
         !---------------------------------------------------------------------!
         character(len=200) :: err_message
         integer(int32) :: io_status
         !---------------------------------------------------------------------!
         open(partial_file_unit, file=trim(partial_xs_file_name),              &
            form='formatted', status='unknown', iostat = io_status,            &
            iomsg = err_message)
         call file_io_status(io_status, err_message, partial_file_unit, "o")
         !---------------------------------------------------------------------!
         call write_message( "  jtot  iblock  v1_f  j1_f  <-  v1_i  j1_i'" //  &
            repeat(" ", 13) // "K.E." // repeat(" ", 16) // "XS",              &
            unit_ = partial_file_unit)
         !---------------------------------------------------------------------!
      end subroutine save_partial_xs_file_header
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine save_partial_xs_single_block(jtot_, block_number_,            &
         open_basis_levels_, xs_block)
         !! ...
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: jtot_
         integer(int32), intent(in) :: block_number_
         integer(int32), intent(in) :: open_basis_levels_(:)
         real(dp), intent(in) :: xs_block(:)
         !---------------------------------------------------------------------!
         character(len=200) :: partial_line
         integer(int32) :: number_of_open_basis_levels_, index_1_, index_2_,   &
            cross_section_index_
         !---------------------------------------------------------------------!
         number_of_open_basis_levels_ = size(open_basis_levels_)
         do index_1_ = 1, number_of_open_basis_levels_
            do index_2_ = 1, number_of_open_basis_levels_
               cross_section_index_ = (index_1_-1)                             &
                  * number_of_open_basis_levels_ + index_2_
               write(partial_line,                                             &
                  "(I6,2X,I6,2X,I4,2X,I4,6X,I4,2X,I4,2X,E16.8,2X,E16.8)")      &
                  jtot_, block_number_,                                        &
                  vib_levels(open_basis_levels_(index_1_)),                    &
                  rot_levels(open_basis_levels_(index_1_)),                    &
                  vib_levels(open_basis_levels_(index_1_)),                    &
                  rot_levels(open_basis_levels_(index_1_)),                    &
                  (total_energy()-internal_energies(open_basis_levels_(index_1_)))&
                     * hartree_to_cm, xs_block(cross_section_index_)
               call write_message(partial_line, unit_ = 12)
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine save_partial_xs_single_block
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      !                           Threshold checking
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine determine_largest_cross_sections(open_basis_levels_,          &
         cross_sections_, max_elastic_cross_section_,                          &
         max_inelastic_cross_section_, max_elastic_index_,                     &
         max_inelastic_index_1_, max_inelastic_index_2_)
         !! Determine the largest partial elastic and inelastic cross-sections
         !! in a given set of cross-sections.
         !---------------------------------------------------------------------!
         integer(int32), intent(in) :: open_basis_levels_(:)
            !! array holding indices to open basis levels
         real(dp), intent(in) :: cross_sections_(:)
            !! array holding partial cross-sections
         real(dp), intent(out) :: max_elastic_cross_section_
            !! largest elastic cross-section
         real(dp), intent(out) :: max_inelastic_cross_section_
            !! largest inelastic cross-section
         integer(int32), intent(out) :: max_elastic_index_
            !! index pointing to the largest elastic cross-section
         integer(int32), intent(out) :: max_inelastic_index_1_
            !! first index pointing to the largest inelastic cross-section
         integer(int32), intent(out) :: max_inelastic_index_2_
            !! second index pointing to the largest inelastic cross-section
         !---------------------------------------------------------------------!
         integer(int32) :: number_of_open_basis_levels_, index_1_, index_2_,   &
            cross_section_index_
         !---------------------------------------------------------------------!
         number_of_open_basis_levels_ = size(open_basis_levels_)
         !---------------------------------------------------------------------!
         ! Initialize output values
         !---------------------------------------------------------------------!
         max_elastic_cross_section_ = 0.0_dp
         max_inelastic_cross_section_ = 0.0_dp
         max_elastic_index_ = 0
         max_inelastic_index_1_ = 0
         max_inelastic_index_2_ = 0
         !---------------------------------------------------------------------!
         ! Iterate over all combinations of open basis levels
         !---------------------------------------------------------------------!
         do index_1_ = 1, number_of_open_basis_levels_
            do index_2_ = 1, number_of_open_basis_levels_
               cross_section_index_ = (index_1_-1)                             &
                  * number_of_open_basis_levels_ + index_2_
               if (open_basis_levels_(index_2_) == open_basis_levels_(index_1_)) then
                  !------------------------------------------------------------!
                  ! Elastic cross-section
                  !------------------------------------------------------------!
                  if (cross_sections_(cross_section_index_)                    &
                     > max_elastic_cross_section_) then
                     max_elastic_cross_section_                                &
                        = cross_sections_(cross_section_index_)
                     max_elastic_index_ = index_1_
                  endif
               else
                  !------------------------------------------------------------!
                  ! Inelastic cross-section
                  !------------------------------------------------------------!
                  if (cross_sections_(cross_section_index_)                    &
                     > max_inelastic_cross_section_) then
                     max_inelastic_cross_section_                              &
                        = cross_sections_(cross_section_index_)
                     max_inelastic_index_1_ = index_1_
                     max_inelastic_index_2_ = index_2_
                  endif
               endif
            enddo
         enddo
         !---------------------------------------------------------------------!
      end subroutine determine_largest_cross_sections
      !------------------------------------------------------------------------!
      !------------------------------------------------------------------------!
      subroutine check_cross_section_thresholds(largest_elastic_xs_,           &
         largest_inelastic_xs_, consecutive_elastic_, consecutive_inelastic_,  &
         terminate_)
         !! Checks if the elastic_xs_threshold (threshold for elastic XS)
         !! and inelastic_xs_threshold (threshold for inelastic XS) conditions
         !! are already fulfilled.
         !---------------------------------------------------------------------!
         real(dp), intent(in) :: largest_elastic_xs_
            !! largest elastic XS in the block
         real(dp), intent(in) :: largest_inelastic_xs_
            !! largest inelastic XS in the block
         integer(int32), intent(inout) :: consecutive_elastic_
            !! number of consecutive blocks meeting condition on elastic XS
         integer(int32), intent(inout) :: consecutive_inelastic_
            !! number of consecutive blocks meeting condition on inelastic XS
         logical, intent(inout) :: terminate_
            !! flag to indicate termination of loop based on thresholds
         !---------------------------------------------------------------------!
         logical :: is_elastic_xs_within_threshold, is_inelastic_xs_within_threshold
         !---------------------------------------------------------------------!
         terminate_ = .false.

         is_elastic_xs_within_threshold                                        &
            = (largest_elastic_xs_ <= elastic_xs_threshold)
         is_inelastic_xs_within_threshold                                      &
            = (largest_inelastic_xs_ <= inelastic_xs_threshold)

         if (is_elastic_xs_within_threshold) then
            consecutive_elastic_ = consecutive_elastic_ + 1
         else
            consecutive_elastic_ = 0
         endif
         
         if (is_inelastic_xs_within_threshold) then
            consecutive_inelastic_ = consecutive_inelastic_ + 1
         else
            consecutive_inelastic_ = 0
         endif

         if ((consecutive_elastic_ >= consecutive_blocks_threshold).and.       &
            (consecutive_inelastic_ >= consecutive_blocks_threshold)) then
            terminate_ = .true.
         endif
         !---------------------------------------------------------------------!
      end subroutine check_cross_section_thresholds
   !---------------------------------------------------------------------------!
end module state_to_state_cross_sections_mod
