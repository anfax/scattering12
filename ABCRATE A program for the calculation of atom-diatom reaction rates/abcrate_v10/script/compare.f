!!!*************************************************************
! 文件/File: compare.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: compare.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      program compare

      implicit none

      real*8        temp,r1,r3,diff
      integer       idx,iunit,maxrat,maxtmp,ncol,ieof,ncolt,i,ntemp
      logical       l1d,l3d
      character*2   cdim
      character*5   cvtst,fname
      character*8   ctun

      parameter(maxtmp=100)
      parameter(maxrat=35)

      dimension temp(maxtmp),r1(maxtmp,maxrat,2,2),
     &          r3(maxtmp,maxrat,2,2),diff(maxtmp,maxrat,2),
     &          cvtst(maxrat),ctun(maxrat)

        do 10 i=1,2

          l1d   = .false.
          l3d   = .false.
          idx   = 1
          ncol  = 0
          ncolt = 0
          ntemp = 0
          iunit = 7
          cdim  = '  '

          if (i .eq. 1) fname = 'file1'
          if (i .eq. 2) fname = 'file2'

          open(unit=iunit,file=fname)

   20     call read1 (iunit,cdim,maxrat,cvtst,ctun,idx,ncol,ntemp,ieof)

            if (ieof .ne. 0) goto 30

            if (cdim .eq. '1D') then
              l1d = .true.
              call read2 (iunit,maxtmp,maxrat,temp,r1,ntemp,idx,ncol,
     &                    i,1,ieof)
            else if (cdim .eq. '3D') then
              l3d = .true.
              call read2 (iunit,maxtmp,maxrat,temp,r3,ntemp,idx,ncol,
     &                    i,2,ieof)
            else
              write(6,*) 'Error - unable to determine calculation type.'
              stop
            end if

            idx=idx+ncol
            if (idx .gt. ncolt) ncolt = idx

          if (ieof .eq. 0) goto 20

   30     continue

          close(unit=iunit)

   10   continue

        ncolt = ncolt - 1

        if (l1d) then
          call cdiff (maxtmp,maxrat,r1,diff,ntemp,ncolt,1)
          call report (maxtmp,maxrat,diff,ntemp,ncolt,1,cvtst,ctun)
        end if

        if (l3d) then
          call cdiff (maxtmp,maxrat,r3,diff,ntemp,ncolt,2)
          call report (maxtmp,maxrat,diff,ntemp,ncolt,2,cvtst,ctun)
        end if

      stop

      end
c
      subroutine read1 (iunit,cdim,maxrat,cvtst,ctun,idx,ncol,ntemp,
     &                  ieof)

      implicit none

      integer       iunit,ieof,idx,ncol,icol,maxrat,i,ntemp
      logical       done
      character*2   cdim
      character*5   cvtst(maxrat)
      character*8   ctun(maxrat)
      character*22  cmp
      character*132 cline

        ieof = 0
        done = .false.
        cmp  = 'Forward rate constants'

   10   read(unit=iunit,fmt=2000,end=500) cline

          if (cline(9:30) .eq. cmp) then
            cdim  = cline(36:37)
            idx   = 1
            ncol  = 0
            ntemp = 0
          end if

          if (cline(10:13) .eq. 'T(K)')  then
            if (ncol .eq. 0) then
              do 20 i=19,91,12
   20           if (cline(i:i+3) .ne. '    ') ncol=(i-7)/12
            end if
            do 30 i=1,ncol
              icol=12*i+7
   30         cvtst(idx-1+i)=cline(icol:icol+4)
            read(unit=iunit,fmt=2000,end=500) cline
            do 40 i=1,ncol
              icol=12*i+7
   40         ctun(idx-1+i)=cline(icol:icol+7)
            done = .true.
          end if

        if (.not. done) goto 10

        goto 900

  500   ieof=1

  900 return

 2000 format(132a)

      end
c
      subroutine read2 (iunit,maxtmp,maxrat,temp,r,ntemp,idx,ncol,
     &                  ifil,idim,ieof)

      implicit none

      integer       iunit,maxtmp,maxrat,ntemp,idx,ncol,ifil,idim,ieof,
     &              i,j
      real*8        tmp,temp,r
      character*132 cline

      dimension  temp(maxtmp),r(maxtmp,maxrat,2,2)

        ieof = 0

        if (ntemp .eq. 0) then
          i=0
   10     read(unit=iunit,fmt=2000,end=500) cline
            if (cline(9:16) .ne. '        ') then
              i=i+1
              read(cline,2010) temp(i),(r(i,idx-1+j,ifil,idim),j=1,ncol)
              goto 10
            end if
        else
          do 20 i=1,ntemp
            read(unit=iunit,fmt=2010,end=500) tmp,
     &           (r(i,idx-1+j,ifil,idim),j=1,ncol)
   20       if (tmp .ne. temp(i)) goto 600
        end if

        if (ntemp .eq. 0) ntemp=i

        goto 900

  500   ieof = 1

        goto 900

  600   write(6,*) 'Input error.  Temperatures do not match in read2.'
        stop

  900   continue

      return

 2000 format(132a)
 2010 format(6x,f10.3,7e12.4)

      end
c
      subroutine cdiff (maxtmp,maxrat,r,diff,ntemp,ncolt,idim)

      real*8  r,diff
      integer ntemp,ncolt,idim,maxtmp,maxrat

      dimension r(maxtmp,maxrat,2,2),diff(maxtmp,maxrat,2)

        do 10 i=1,ntemp
          do 10 j=1,ncolt
   10       diff(i,j,idim)=(r(i,j,2,idim)-r(i,j,1,idim))/r(i,j,1,idim)

      return

      end
c
      subroutine report (maxtmp,maxrat,diff,ntemp,ncolt,idim,cvtst,ctun)

      implicit none

      real*8       diff,sum,rmserr
      integer      maxtmp,maxrat,ntemp,ncolt,idim,i,j,ilen,iloc,strlen
      character*5  cvtst
      character*8  ctun
      character*14 cout

      dimension diff(maxtmp,maxrat,2),cvtst(maxrat),ctun(maxrat)

        if (idim .eq. 1) then
          write(6,1000) '1D'
        else
          write(6,1000) '3D'
        end if

        write(6,1010)

        do 10 i=1,ncolt
          sum=0.d0
          do 20 j=1,ntemp
   20       sum=sum+diff(j,i,idim)*diff(j,i,idim)
          rmserr=dsqrt(sum/dble(ntemp))
          cout='            '
          ilen=strlen(cvtst(i),5)
          cout(1:ilen)=cvtst(i)(1:ilen)
          iloc=ilen+1
          ilen=strlen(ctun(i),8)
          if (ilen .gt. 0) then
            cout(iloc:iloc)='/'
            cout(iloc+1:iloc+ilen)=ctun(i)(1:ilen)
          end if
          write(6,1020) cout,rmserr
   10   continue

      return

 1000 format(/,1x,'Results for ',a2,' reaction',/)
 1010 format(1x,'Method',14x,'RMS Error (relative)')
 1020 format(1x,a14,5x,1p,d11.4)

      end
c
      integer function strlen (str,idlen)

      implicit none

      integer     idlen,itmp
      character*1 str

      dimension   str(idlen)

        itmp=idlen
   10   if (str(itmp) .ne. ' ') goto 20
          itmp=itmp-1
          if (itmp .le. 0) goto 20
        goto 10

   20   strlen=itmp

      return

      end
