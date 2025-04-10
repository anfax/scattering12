!!!*************************************************************
! 文件/File: main_potini.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: main_potini.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

      program potini
!     *********************************************************************
!     Initialize the potential, reactant and product wf, and dipole transition moment
!     to be read to proceed for several wave packet calculations witj MadWave3
!     for AB+C --> AC + C
!     for 01+2 --> 02 + 1
!     for different initial states
!     *********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_Hphi_01y2
      
      implicit none
      character*40 filename
      integer ierr,ielec

      include "mpif.h"

      
c! Initialize MPI environment and get proc's ID and number of proc in
c! the partition.

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      write(filename,'("salpot."i3.3)')idproc
      open(6,file=filename,status='unknown')
      write(6,'(40("="),//)')
      write(6,'(10x,"Potini for Madwave3 version 6 ",//)')
      write(6,*)' output of proc. idproc= ',idproc,' of nproc= ',nproc
      write(6,'(/,40("="),//)')
      
!     initialization of data


      call input_grid
      call pot0
      call  paralelizacion
      write(6,*)
      write(6,*)'  --- Initialization of the potential ---'
      write(6,*)

      nelec=nelecmax
      call setxbcpotele(iomdiat,iomatom,sigdiat,sigatom
     &                     ,nelec,nelecmax)
      write(6,*)
      write(6,*)'  --- end initialization of the potential ---'
      write(6,*)
      if(nelec.gt.nelecmax.or.nelec.eq.0)then
          write(6,*)'  !!! nelec= ',nelec
     &                            ,'  while nelecmax= ',nelecmax
          call flush(6)
          stop
      endif
         open(10,file='pot/cont.pot',status='new')
         write(10,*)nelec
         do ielec=1,nelec
            write(10,*)iomdiat(ielec),iomatom(ielec)
     &           ,sigdiat(ielec),sigatom(ielec)
         enddo
         close(10)

!     determining basis

      call basis
      
! reactants functions calculation

      call angular_functions

      if(npun1.eq.1)then

         write(6,*)' --> calling rigid rotor energies <-- '
         call flush(6)
         call radial_functions01eq_write

      else
         call radial_functions01_write

         if(iprod.eq.2)then
 
            call product_radialf_write

         endif
      endif

!     determining potential

      call pot1

!     determining electric dipole transition

      if(iphoto.gt.0.and.iphoto.le.2)then
         call pot2                ! to initialize global indexes
         if(iphoto.eq.1)then
            call set_trans_dipole
         endif
      elseif(iphoto.eq.3)then
         call pot2                ! to initialize global indexes
         call set_coupling        ! for electronic predissociation
      endif
      
      stop
      end program
