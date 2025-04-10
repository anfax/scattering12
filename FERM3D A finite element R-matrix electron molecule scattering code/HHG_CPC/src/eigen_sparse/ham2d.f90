!!!*************************************************************
! 文件/File: ham2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: ham2d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************



SUBROUTINE ham2d(h2d,h1,o1,h2,o2,nb1,nb2)
  implicit none

  INTEGER :: nb1, nb2
  REAL(kind=8) :: h2d(-60:60,*)
  REAL(kind=8) :: h1(-5:5,nb1), o1(-5:5,nb1)
  REAL(kind=8) :: h2(-5:5,nb2), o2(-5:5,nb2)

  INTEGER :: col2, d2, col1, dstrt, cstrt, md1, md2, strd1, strd2
  REAL(kind=8) :: tol, Otmp, Htmp

  tol = 1.d-20
  md1 = min(5,nb1-1)
  md2 = min(5,nb2-1)
  strd1 = 2*md1 + 1
  strd2 = MIN(nb1,strd1)

  DO col2 = 1,nb2
    cstrt = (col2-1)*nb1
    DO d2 = -md2,md2
      Otmp = o2(d2,col2)
      Htmp = h2(d2,col2)
      dstrt = d2*strd2
      IF (abs(Otmp) > tol) THEN
        DO col1 = 1,nb1
          h2d(dstrt-md1:dstrt+md1,cstrt+col1) = h2d(dstrt-md1:dstrt+md1,cstrt+col1) &
          &  + Otmp*h1(-md1:md1,col1)
        END DO
      ENDIF
      IF (abs(Htmp) > tol) THEN
        DO col1 = 1,nb1
          h2d(dstrt-md1:dstrt+md1,cstrt+col1) = h2d(dstrt-md1:dstrt+md1,cstrt+col1) &
          &  + Htmp*o1(-md1:md1,col1)
        END DO
      ENDIF
    END DO
  END DO

  RETURN
END SUBROUTINE ham2d
