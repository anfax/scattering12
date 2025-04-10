!!!*************************************************************
! 文件/File: sparsemultiply.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: sparsemultiply.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


SUBROUTINE sparsemultiply(product,vector,matvals,rowind,colptr,ntot,nnz,nvec)
  implicit none

  INTEGER, intent(IN) :: ntot,nvec,nnz
  INTEGER, intent(IN) :: rowind(nnz), colptr(0:ntot)
  REAL(kind=8), intent(OUT) :: product(ntot*nvec)
  REAL(kind=8), intent(IN)  :: vector(ntot*nvec), matvals(nnz)
  INTEGER :: col,j,row,strt,iv,shft

  product(:) = 0.d0

  shft = 1 - colptr(0)
  strt = 0
  DO iv = 1,nvec
    DO col = 1,ntot
      DO j = colptr(col-1)+shft,colptr(col)+shft-1
        row = rowind(j)+shft
        product(strt+row) = product(strt+row) + matvals(j)*vector(strt+col)
      END DO !j
    END DO !col
    strt = strt + ntot
  END DO

  RETURN
END SUBROUTINE sparsemultiply

