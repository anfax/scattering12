!!!*************************************************************
! 文件/File: plgndr.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: plgndr.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

module leg
implicit none

contains
 
FUNCTION plgndr(l,m,x)
      INTEGER*8 l,m,s
      REAL*8 plgndr,x
      INTEGER*8 i,ll
      REAL*8 fact1,pll,pmm,pmmp1,somx2
      s = 0
     if(m.lt.0) then
  m=abs(m)
  s=1
     endif
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact1=1.
        do 11 i=1,m
          pmm=-pmm*fact1*somx2
          fact1=fact1+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      if(s.eq.1)then
  plgndr=((-1)**m)*(real(fact(l-m))/fact(l+m))*plgndr
  m=-m
    end if
      return
      END FUNCTION plgndr


end program triatomic

end module
