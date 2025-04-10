!!!*************************************************************
! 文件/File: SuperLU.bridge.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: SuperLU.bridge.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

module pardiso_mod
use nrtype, only : i4b,dbl
contains

subroutine pardiso_sym(Gamma0sparse,rowind,Colptr,nnz,nClosed,nOpen,RHS_co,ierr,mtype)
  integer(kind=i4b),intent(in) :: nnz,nClosed,nOpen,mtype
  integer(kind=i4b),intent(inout):: ierr
  INTEGER(kind=i4b) :: maxfct, mnum,  phase, error, msglvl,LDB
  INTEGER(kind=i4b),intent(inout) :: Colptr(nClosed+1), rowind(nnz)
  real(kind=dbl),intent(inout) :: Gamma0sparse(nnz),RHS_co(nClosed,nOpen)
LDB=NClosed
write(6,*)'before SuperLU'
 call dslusolve(Nclosed,nnz,Nopen,Gamma0sparse,rowind,colptr,RHS_co,LDB,ierr,'1',3)
write(6,*)'after SuperLU'

end subroutine pardiso_sym
end module pardiso_mod
