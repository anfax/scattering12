!!!*************************************************************
! 文件/File: multipl_norm.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: multipl_norm.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


SUBROUTINE multipl_norm(product,vector,matvals,rowind,ntot,nvec,dime)
  implicit none

  INTEGER, intent(IN) :: ntot,nvec,dime
  INTEGER, intent(IN) :: rowind(ntot,dime)
  REAL(kind=8), intent(OUT) :: product(nvec,nvec)
  REAL(kind=8), intent(IN)  :: vector(ntot,nvec), matvals(ntot,dime)
  real(kind=8),dimension(:,:),allocatable :: product_interm
  INTEGER :: col,j,row,strt,iv,shft,i,k,l
  allocate(product_interm(ntot,nvec))
  product(:,:) = 0.d0
  product_interm(:,:)=0.d0
write(6,*)'inside multipl_norm'
  DO i = 1,ntot,8
    DO k = 1,nvec
      do l=1,ntot,8
        DO j = 1,dime
          if (rowind(i,j).eq.l) then
            product_interm(i,k) = product_interm(i,k) + matvals(i,j)*vector(l,k)
            exit
          end if
        end do
      END DO !j
    END DO !col
  END DO
write(6,*)'after 1st product'
!  product= - matmul(transpose(product_interm),vector)

! Call Blas 3 subroutine to perform this multiplication
call dgemm ( 'T', 'N', nvec, nvec, ntot, 1.d0, product_interm, ntot, vector, ntot, &
     &                   0.d0, product, nvec )
write(6,*)'after 2nd product'
deallocate(product_interm)
  RETURN
END SUBROUTINE multipl_norm


subroutine Normalize(ntot,nvec,cc)
implicit none
  real(kind=8),dimension(:,:),allocatable :: Overlap, rowind,product
  real(kind=8),dimension(ntot,nvec) :: cc
  integer(kind=4) :: ntot,nvec,dime=216,ind,ind_max,i,j,colind
  allocate(Overlap(ntot,dime))
  allocate(rowind(ntot,dime))
  allocate(product(nvec,nvec))
        !*COMM* Restore Part
rowind(:,:)=0.d0
Overlap(:,:)=0.d0
product(:,:)=0.d0
        read(1000,*)ind_max
        write(6,*)'before read'
        do i=1,ind_max
        read(1000,*)ind,colind
           do j=1,colind
                read(1000,*)Overlap(i,j),rowind(i,j)
           end do
        end do
call multipl_norm(product,cc,Overlap,rowind,ntot,nvec,dime)
write(6,*)'after multipl_norm'
  deallocate(Overlap)
  deallocate(rowind)
  deallocate(product)

end subroutine Normalize
