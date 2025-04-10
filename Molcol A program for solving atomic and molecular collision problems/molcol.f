!!!*************************************************************
! 文件/File: molcol.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: molcol.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

      program molcol
      implicit double precision (a-h,o-z)
c
c     This program is based on the following publications:
c     Launay J.-M. 1977 j. phys. b: atom. molec. phys. 10, 3665
c     Flower D. R. and Launay J.-M. 1977 j. phys. b: atom. molec. phys.
c                                        10, 3673
c     Flower D. R. and Launay J.-M. 1977 j. phys. b: atom. molec. phys.
c                                        10, L229
c     Launay J.-M. 1976 j. phys. b: atom. molec. phys. 9, 1823
c     please cite Launay (1977), in addition to computer physics
c     communications, in any publication involving the use of this
c     program, called "molcol"
c
c     dimension super-array x for dynamic core allocation
c
      parameter (nmmx=4200000)
      dimension x(nmmx),n(64)
c
c     open the input and output files : dat1 specifies dimensions and
c     energies ; dat2 contains the potential
c
      open (5, file = 'molcol.dat1', status = 'unknown')
      open (9, file = 'molcol.dat2', status = 'unknown')
      open (6, file = 'molcol.out', status = 'unknown')
c
c     maximum total size of dimensioned variables; dimensions allocated
c     in dime, where actual total size is determined and printed
c     present example allows for a total of 4 200 000 real*8
c     array elements, equivalent to a memory requirement of 34 mbytes
c
      call dime (n,nmmx)
c
c     call prncp, dynamically allocating the storage
c
      call prncp (x(n(1 )),x(n(2 )),x(n(3 )),x(n(4 )),x(n(5 )),x(n(6 ))
     &           ,x(n(7 )),x(n(8 )),x(n(9 )),x(n(10)),x(n(11)),x(n(12))
     &           ,x(n(13)),x(n(14)),x(n(15)),x(n(16)),x(n(17)),x(n(18))
     &           ,x(n(19)),x(n(20)),x(n(21)),x(n(22)),x(n(23)),x(n(24))
     &           ,x(n(25)),x(n(26)),x(n(27)),x(n(28)),x(n(29)),x(n(30))
     &           ,x(n(31)),x(n(32)),x(n(33)),x(n(34)),x(n(35)),x(n(36))
     &           ,x(n(37)),x(n(38)),x(n(39)),x(n(40)),x(n(41)),x(n(42))
     &           ,x(n(43)),x(n(44)),x(n(45)),x(n(46)),x(n(47)),x(n(48))
     &           ,x(n(49)),x(n(50)),x(n(51)),x(n(52)),x(n(53)),x(n(54))
     &           ,x(n(55)),x(n(56)),x(n(57)),x(n(58)),x(n(59)),x(n(60))
     &           ,x(n(61)),x(n(62)),x(n(63)),x(n(64)))
      stop
      end
c-----------------------------------------------------------------------
      subroutine dime (n,nmmx)
c
c     called from the main program
c     determines the dimensions of the principal arrays for dynamic
c     core allocation
c
      implicit double precision (a-h,o-z)
c
c     This routine computes  array dimensions
c     nptmx  :  maximum number of potential curve points
c     n1mx   :  maximum number of (j1) states for particle 1
c     n2mx   :  maximum number of (j2) states for particle 2
c     nvmx   :  maximum number of channels
c     nj12mx :  maximum number of j12 values when coupling j1 and j2
c     ncjmx  :  maximum number of elements in BF-SF unitary transformation
c     nglmmx :  maximum of algebraic coefficients
c     nq12mx :  maximum number of potential curves
c
      common /ndim/ nptmx,nvmx,njsmx,nj12mx,nq12mx,ndgnmx,ncjmx,nglmmx
     &             ,nblkmx,n1mx,n2mx
      dimension n(*)
      namelist /dmnsns/ n1mx,n2mx,nvmx,nj12mx,ncjmx,nglmmx,nptmx,nq12mx
      nblkmx = 100
c
c     dimensions read from data file molcol.dat1 and written to output file
c     molcol.out
c     subsequent checks lead to the printing of error messages if nvmx,
c     ncjmx, nglmmx, nptmx, nq12mx, nblkmx, or the overall dimension ndmmx
c     is exceeded
c
      read (5,dmnsns)
      write (6, 9000)
      write (6,dmnsns)
      njsmx = n1mx * n2mx
      njsm2 = njsmx / 2 + 1
      nj12m2 = nj12mx / 2 + 1
      nvmx2 = nvmx * nvmx
      nvmxs2 = nvmx / 2 + 1
c
      n(1)  = 1
      n(2)  = n(1)  + nj12mx
      n(3)  = n(2)  + nj12mx
      n(4)  = n(3)  + nj12mx
      n(5)  = n(4)  + nj12mx
      n(6)  = n(5)  + nj12mx
      n(7)  = n(6)  + nj12mx
      n(8)  = n(7)  + nj12mx
      n(9)  = n(8)  + nq12mx
      n(10) = n(9)  + nq12mx
      n(11) = n(10) + nq12mx
      n(12) = n(11) + nq12mx * nj12mx * (nj12mx+1) / 2
      n(13) = n(12) + nj12mx
      n(14) = n(13) + nj12mx
      n(15) = n(14) + nvmx
      n(16) = n(15) + nvmx
      n(17) = n(16) + nvmx
      n(18) = n(17) + nvmxs2
      n(19) = n(18) + nvmxs2
      n(20) = n(19) + nq12mx
      n(21) = n(20) + nq12mx
      n(22) = n(21) + nq12mx
      n(23) = n(22) + nvmx2
      n(24) = n(23) + nvmx2
      n(25) = n(24) + nvmx2
      n(26) = n(25) + nvmx2
      n(27) = n(26) + nvmx2
      n(28) = n(27) + nvmx
      n(29) = n(28) + nvmx
      n(30) = n(29) + nvmx
      n(31) = n(30) + nvmx
      n(32) = n(31) + nvmx
      n(33) = n(32) + nvmx
      n(34) = n(33) + nvmx
      n(35) = n(34) + nptmx
      n(36) = n(35) + nptmx
      n(37) = n(36) + nptmx
      n(38) = n(37) + nq12mx
      n(39) = n(38) + nq12mx
      n(40) = n(39) + njsmx
      n(41) = n(40) + njsmx
      n(42) = n(41) + njsmx
      n(43) = n(42) + njsmx
      n(44) = n(43) + njsmx
      n(45) = n(44) + njsmx
      n(46) = n(45) + njsmx
      n(47) = n(46) + nvmxs2
      n(48) = n(47) + nvmxs2
      n(49) = n(48) + nvmxs2
      n(50) = n(49) + nglmmx
      n(51) = n(50) + nvmxs2
      n(52) = n(51) + nvmxs2
      n(53) = n(52) + ncjmx
      n(54) = n(53) + nvmx + nvmx
      n(55) = n(54) + nvmxs2
      n(56) = n(55) + nvmxs2
      n(57) = n(56) + nvmxs2
      n(58) = n(57) + nq12mx
      n(59) = n(58) + nj12mx * nj12mx
      n(60) = n(59) + njsm2
      n(61) = n(60) + njsm2
      n(62) = n(61) + njsmx
      n(63) = n(62) + njsmx * njsmx
      n(64) = n(63) + njsmx * njsmx
      nfin  = n(64) + nvmx2
c
c     print the total the number of storage locations used nfin and
c     the maximum number that can be allocated nmmx
c
      write (6, 9010) nfin,nmmx
c
c     program stops if insufficient storage allocated in main
c
      if (nfin .gt. nmmx) stop
      return
 9000 format (' *****   DIMENSIONS  USED  *****')
 9010 format ('   REAL*8 WORDS USED : ',i8
     &       ,' ; REAL*8 WORDS AVAILABLE : ',i8)
      end
c-----------------------------------------------------------------------
      subroutine prncp (fl1jj ,fs1jj ,fj1jj ,fl2jj ,fs2jj ,fj2jj
     &                 ,fj12  ,q1q   ,q2q   ,q12   ,xred  ,pprt
     &                 ,ekjj  ,xmg   ,fl    ,xnorm ,idxbf ,idxsf
     &                 ,q1m   ,q2m   ,fmu   ,aa    ,bb    ,cc
     &                 ,dd    ,ee    ,enbf  ,ak    ,akr   ,sb
     &                 ,sbp   ,sn    ,snp   ,rpot  ,vpot  ,cm
     &                 ,nptpot,ndbpot,fl1js ,fs1js ,fj1js ,fl2js
     &                 ,fs2js ,fj2js ,fdgn  ,nsom  ,ndbsom,ivv
     &                 ,glm   ,ndmblk,ndbmat,cj    ,cr    ,ir1p
     &                 ,ir2p  ,idbglm,vv    ,vmat  ,idx1rg,idx2rg
     &                 ,dsgma ,opc   ,sct   ,zz)
      implicit double precision (a-h,o-z)
c
c     variables relating to the numerical integration of the coupled
c     equations transferred via common block /intgr/
c
      common /intgr/ rdb,rfn,fpt,enfdml,rmin,stbmx,fllmx,enmin,nstbmx
      common /steps/ step,eneff,enstb,deltae
c
c     integer variables which determine the level of printout transferred
c     via common block /ecrit/
c
      common /ecrit/ ikmat,it2,iopc,isctp,isctot,itest,iradin
c
c     maximum dimensions transferred via common block /ndim/
c
      common /ndim/ nptmx,nvmx,njsmx,nj12mx,nq12mx,ndgnmx,ncjmx,nglmmx
     &             ,nblkmx,n1mx,n2mx
c
c     data relating to the dispersion energy transferred via common/disp/
c
c      common /disp/ xdsp(200),cdsp(200),cmdsp(200),nptdsp(25)
c     &             ,ndbdsp(25)
c     reduced mass of projectile and target rdmas
c
      common /mas/ rdmas
c
c     data relating to vibrational coupling transferred via common/vib/
c     the following limits are implied by the following dimensions:
c     nq12mx.le.150, nj12mx.le.90, nglmmx.le.400000
c
      common /vib/ fn2q(150),fnp2q(150),fn1jj(90),fn2jj(90),fj2jjj(90)
c
c     2 total parities par, 20 collision energies en (in k),
c     30 projectile e1k and target e2k energy levels (implies n1mx,
c     n2mx.le.30 and n1mx*n2mx.le.900)
c     levels identified by their energies e1k, e2k (in k), their vibrational
c     quantum numbers fn1j, fn2j and the values of l, s and j (fl1j, fl2j;
c     fs1j, fs2j; fj1, fj2), where j = l + s (vector sum)
c     note the conversion factor: 1 cm^(-1) = 1.4388 k = 1.2399e-4 ev
c
      dimension par(2),en(20), e1k(30),e2k(30), fn1j(30),fl1j(30)
     &         ,fs1j(30), fj1(30), fn2j(30),fl2j(30),fs2j(30),fj2(30)
     &         ,fn1js(900),fn2js(900)
c
c      nvmx.le.801 is implied by the following dimension statement:
c
      dimension trp(801),trt(801),sctinp(801),sctint(801)
     &         ,ipiv(801),work(801)
      dimension
     &     fl1jj (*),fs1jj (*),fj1jj (*),fl2jj (*),fs2jj (*),fj2jj (*)
     &    ,fj12  (*),q1q   (*),q2q   (*),q12   (*),xred  (*),pprt  (*)
     &    ,ekjj  (*),xmg   (*),fl    (*),xnorm (*),idxbf (*),idxsf (*)
     &    ,q1m   (*),q2m   (*),fmu   (*),aa    (*),bb    (*),cc    (*)
     &    ,dd    (*),ee    (*),enbf  (*),ak    (*),akr   (*),sb    (*)
     &    ,sbp   (*),sn    (*),snp   (*),rpot  (*),vpot  (*),cm    (*)
     &    ,nptpot(*),ndbpot(*),fl1js (*),fs1js (*),fj1js (*),fl2js (*)
     &    ,fs2js (*),fj2js (*),fdgn  (*),nsom  (*),ndbsom(*),ivv   (*)
     &    ,glm   (*),ndmblk(*),ndbmat(*),cj    (*),cr    (*),ir1p  (*)
     &    ,ir2p  (*),idbglm(*),vv    (*),vmat  (*),idx1rg(*),idx2rg(*)
     &    ,dsgma (*),opc   (*),sct   (*),zz    (*)
      namelist /entree/ nj1,fn1j,fl1j,fs1j,fj1,nj2,fn2j,fl2j,fs2j,fj2
     &                 ,e1k,e2k,fmas1,fmas2
c
c     nj1   : nombre of (j1) sates for particule 1
c     fn1j  : values of the supplementary quantum number
c     fl1j  : values of l1
c     fs1j  : values of s1
c     fj1   : values of j1
c     e1k   : energies of (j1) sates in degres kelvin
c     fmas1 : masse of particule 1 ( proton mass = 1)
c     fmas2 : masse of particule 2
c
      namelist /pot/ nq12,nnq12
c
c     nnq12 : number of potential curves
c     nq12  : number of values of q12 = q1 + q2 (vector sum)
c
      namelist /intg/ fjdb,fjfn,fjps,en,nenrg,rdb,rfn,fpt,opcmn,rmin
     &               ,deltae,nstbmx,stbmx,integrator
c
c     fjdb   : starting value of total angular momentum
c     fjfn   : last     value of total angular momentum
c     fjps   : increment for total angular momentum
c     nenrg  : nombre of total energies
c     en     : total energies (degrees kelvin)
c     rdb    : first integration point (in a0)
c     rfn    : last  integration point (in a0)
c     fpt    : nomber of points per arch
c              recommended value : between 5 et 10
c     opcmn  : stopping criterion on maximal value of inelastic probas.
c     rmin   : position of minimum of most attractive potential curve
c     deltae : value of minimum (u.a.)
c     nstbmx : nombre maximal de points sans stabilisation
c     stbmx  : stabilisation criterion (between 5 and 10)
c     integrator : .le. 0, de vogelaere method
c                  .eq. 1, johnson method
c                  .eq. 2, johnson-manolopoulos method
c                  .eq. 3, johnson method
c                           with symmetric matrix inverter
c                  .ge. 4, johnson-manolopoulos method
c                           with symmetric matrix inverter
c
      namelist /ecr/ ikmat,it2,iopc,isctp,isctot,itest,iradin,itr,ipot
     &              ,ialg
c
c     indices for printing (0 : no print)
c     itest prints integration tests
c     itr : prints total inelasticities and energy tranfers
c     iradin = 0 :  minimal information
c
      data fmprt,cnvkua /1836.104d0,3.16647d-6/
      data pi /3.141592654d0/
      data par /-1.d0,+1.d0/
      data eps /0.01d0/
      xl(xxx) = 2.*xxx+1.
c
c     check if dimensions exceeded, write error message and stop if
c     it is the case
c
      if (nq12mx .gt. 150) then
         write (6, 9020) nq12mx
         stop
      else if (nj12mx .gt. 90) then
         write (6, 9030) nj12mx
         stop
      else if (n1mx .gt. 30) then
         write (6, 9050) n1mx
         stop
      else if (n2mx .gt. 30) then
         write (6, 9060) n2mx
         stop
      else if (nvmx .gt. 801) then
         write (6, 9070) nvmx
         stop
      endif
c
c     input data read from molcol.dat1
c
      read (5,entree)
      rdmas = (fmas1 * fmas2 / (fmas1+fmas2)) * fmprt
      read (5,pot)
c
c     check if dimension exceeded
c
      if (nnq12 .gt. nq12mx) then
         write (6, 9000) nnq12,nq12mx
         stop
      endif
      read (5,ecr)
      if (iradin .le. 0) then
         ikmat = 0
         it2 = 0
         iopc = 0
         itr = 0
         isctp = 0
         isctot = 0
         itest = 0
      endif
c
c     write (6,1270) nnq12
c     quantum numbers for potential energy curves read from molcol.dat2
c
      do i = 1,nnq12
	 read (9,9300) fn2q(i),fnp2q(i), q1q(i),q2q(i),q12(i)
     &                ,q1m(i),q2m(i),fmu(i)
c
c     fn2q, fnp2q  : supplementary quantum numbers
c     q1m,q2m,fmu : quantum numbers for BF expansion
c     q1q,q2q,q12 : quantum numbers for SF expansion
c
      enddo
      read (5,intg)
c
c     print the input data for confirmation
c
      do i = 1, nj1
         write (6, 9310) fn1j(i),fl1j(i),fs1j(i),fj1(i)
      enddo
      do i = 1, nj2
         write (6, 9310) fn2j(i),fl2j(i),fs2j(i),fj2(i)
      enddo
      write (6, 9310) fmas1,fmas2
      write (6,intg)
c
c     levels are sorted by increasing energy
c
      nj1nj2 = nj1 * nj2
      ekmn = 1.d72
      ekmx  = -1.d72
      do i1 = 1,nj1
         do i2 = 1,nj2
            ek = e1k(i1) + e2k(i2)
            if (ek .le. ekmn) ekmn = ek
            if (ek .ge. ekmx) ekmx = ek
	 enddo
      enddo
      i = 0
 100  do i1 = 1,nj1
         do i2 = 1,nj2
            ek = e1k(i1) + e2k(i2)
            if (ek .eq. ekmn) then
               i = i + 1
               idx1rg(i) = i1
               idx2rg(i) = i2
	    endif
	 enddo
      enddo
      if (ekmn .ne. ekmx) then
         ekmn1 = 1.d72
         do i1 = 1, nj1
            do i2 = 1, nj2
               ek  = e1k(i1) + e2k(i2)
               if (ek .le. ekmn1 .and. ek .gt. ekmn)  ekmn1 = ek
	    enddo
	 enddo
         ekmn = ekmn1
	 go to 100
      endif
c
c     vectorial coupling of j1 and j2 gives j12
c
      nj12 = 0
      do i = 1,nj1nj2
         i1 = idx1rg(i)
         i2 = idx2rg(i)
         fn1js(i) = fn1j(i1)
         fl1js(i) = fl1j(i1)
         fs1js(i) = fs1j(i1)
         fj1js(i) = fj1(i1)
         fn2js(i) = fn2j(i2)
         fl2js(i) = fl2j(i2)
         fs2js(i) = fs2j(i2)
         fj2js(i) = fj2(i2)
         fjmn = dabs(fj1(i1)-fj2(i2))
         fjmx = fj1(i1) + fj2(i2)
         fj = fjmn
 150     continue
         nj12 = nj12 + 1
         fn1jj(nj12) = fn1j(i1)
         fn2jj(nj12) = fn2j(i2)
         fl1jj(nj12) = fl1j(i1)
         fl2jj(nj12) = fl2j(i2)
         fs1jj(nj12) = fs1j(i1)
         fs2jj(nj12) = fs2j(i2)
         fj1jj(nj12) = fj1(i1)
         fj2jj(nj12) = fj2(i2)
         fj2jjj(nj12) = fj2(i2)
         fj12(nj12) = fj
         ekjj(nj12) = e1k(i1) + e2k(i2)
         fj = fj + 1
         if (fj .lt. (fjmx+eps)) go to 150
      enddo
      fj12mx = 0.
      do i = 1, nj12
         if (fj12(i) .gt. fj12mx) fj12mx = fj12(i)
      enddo
      do i = 1,nj1nj2
         write (6,9310) fn1js(i),fl1js(i),fs1js(i),fj1js(i)
     &                 ,fn2js(i),fl2js(i),fs2js(i),fj2js(i)
      enddo
      write (6,9320)
      do i = 1, nj12
c
c     print the list of channels, including the intermediate
c     vector sum j12 = j1 + j2 and the energies (in K)
c
         write (6,9310) fn1jj(i),fl1jj(i),fs1jj(i),fj1jj(i)
     &                 ,fn2jj(i),fl2jj(i),fs2jj(i),fj2jj(i)
     &                 ,fj12(i),ekjj(i)
      enddo
      do i = 1,nj12
         nx = fl1jj(i) + fl2jj(i) + eps
	 pprt(i) = 1.d0
         if (mod(nx,2) .gt. 0) pprt(i) = -pprt(i)
      enddo
      write (6, 9350)
      do i = 1, nnq12
c
c     print the quantum numbers which label the potential energy curves
c     q1m,q2m,fmu relate to the BF representation
c     q1q,q2q,q12 relate to the SF representation
c
	 write (6, 9310) fn2q(i),fnp2q(i)
     &                  ,q1q(i),q2q(i),q12(i),q1m(i),q2m(i),fmu(i)
      enddo
c
c     reads the potential (molcol.dat2)
c     distances and energy in atomic units
c
c     potential is written as v(hartree-fock) + v (dispersion)
c     this formal separation is not always necessary or desirable
c     in the present version of the program, the "dispersion"
c     contribution is commented out and the "hartree-fock" contribution
c     is taken equal to the total energy
c     the potential expansion coefficients (body-fixed frame) are input
c     on a grid of values of the intermolecular distance (distance
c     between the molecular centres of mass) for subsequent cubic
c     spline interpolation
c     optional printout of the potential expansion coefficients
c
      ndb = 1
      do iq = 1,nnq12
         ndbpot(iq) = ndb
         read (9,9200) npnt,(rpot(ndb+i-1),vpot(ndb+i-1),i=1,npnt)
         ndb = ndb + npnt
         npt = ndb
         nptpot(iq) = npnt
      enddo
c
c     check if dimension exceeded
c
      if (npt .gt. nptmx) then
         write (6, 9010) npt,nptmx
         stop
      endif
      do iq = 1,nnq12
         ndb = ndbpot(iq)
         npnt = nptpot(iq)
         if (ipot .ge. 1) write (6, 9210)
     &        npnt,(rpot(ndb+i-1),vpot(ndb+i-1),i=1,npnt)
      enddo
c
      i = 0
      do iq = 1,nnq12
         npnt = nptpot(iq)
         if (npnt .ge. 1) i = i + npnt
      enddo
      npot = i
      if (ipot .gt. 0) write (6, 9500) (vpot(i),i=1,npot)
      id = 0
c
c     spline of potentiel
c
      do iq = 1,nnq12
         ndb = ndbpot(iq)
         npnt = nptpot(iq)
         call spline (npnt,rpot(ndb),vpot(ndb),cm(ndb))
      enddo
c
c     calculation of reduced element (voir Launay 1977)
c
      ixred = 0
      do i = 1, nj12
c
c     function xl(x) evaluates 2*x+1
c
         r2 = xl(fl1jj(i)) * xl(fl2jj(i)) * xl(fj1jj(i))
     &      * xl(fj2jj(i)) * xl(fj12(i))
         do ip = 1,i
            r3 = xl(fl1jj(ip)) * xl(fl2jj(ip)) * xl(fj1jj(ip))
     &         * xl(fj2jj(ip)) * xl(fj12(ip))
            do iq = 1,nq12
               ixred = ixred + 1
c
c     (2*q12+1)**1/2 factor is added in reduced element
c
               r1 = xl(q1q(iq)) * xl(q2q(iq)) * xl(q12(iq))
     &            * xl(q12(iq))
	       r4 = f3j (fl1jj(i),q1q(iq),fl1jj(ip),0.d0,0.d0,0.d0)
c
c     3j, racah and 9j coefficients are evaluated
c
	       r5 = f3j (fl2jj(i),q2q(iq),fl2jj(ip)
     &                  ,0.d0,0.d0,0.d0)
	       r6 = f6j (fl1jj(i),fj1jj(i),fs1jj(i)
     &                  ,fj1jj(ip),fl1jj(ip),q1q(iq))
	       r7 = f6j (fl2jj(i),fj2jj(i),fs2jj(i)
     &                  ,fj2jj(ip),fl2jj(ip),q2q(iq))
	       r8 = f9j (fj1jj(i),fj2jj(i),fj12(i)
     &                  ,fj1jj(ip),fj2jj(ip),fj12(ip)
     &                  ,q1q(iq),q2q(iq),q12(iq))
               nx = fs1jj(i) + fj1jj(ip) + fs2jj(i) + fj2jj(ip)
     &            + q1q(iq) + q2q(iq) + eps
	       sgn = 1.d0
               if (mod(nx,2) .gt. 0) sgn = -sgn
               xred(ixred) = sgn * dsqrt(r1*r2*r3) * r4 * r5 * r6 * r7
     &                     * r8
	    enddo
	 enddo
      enddo
      nxred = ixred
      if (ialg .gt. 0) write (6, 9500) (xred(i),i=1,nxred)
c
c     conversion to atomic units and multiplication by 2*reduced mass
c
      deltae = deltae * rdmas  * 2.d0
      dmskua = 2.d0 * rdmas * cnvkua
      do i = 1,nj1
         e1k(i) = e1k(i) * dmskua
      enddo
      do i = 1,nj2
         e2k(i) = e2k(i) * dmskua
      enddo
      do i = 1,nj12
         ekjj(i) = ekjj(i) * dmskua
      enddo
c
c     loop over collision energy
c
      do ien = 1,nenrg
         write (6, 9400) en(ien)
         en(ien) = en(ien) * dmskua
         do ip = 1,njsmx
            do ipp = 1,njsmx
               sct(njsmx*(ip-1)+ipp) = 0
	    enddo
	 enddo
         enfdml = en(ien) - ekjj(1)
c
c     hmx :  step before potential minimum
c
         fjtot = fjdb - fjps
 200     continue
c
c     loop over partial waves
c
         fjtot = fjtot + fjps
c
c     option to print the value of the total angular momentum quantum
c     number ("PARTIAL WAVE J")
c
         if (iradin .gt. 0) write (6, 9410) fjtot
         nit = 0
         do it = 1,nj1nj2
            i1 = idx1rg(it)
            i2 = idx2rg(it)
            xak = en(ien) - e1k(i1) - e2k(i2)
            if (xak .ge. 0) then
               nit = nit + 1
c
c     cross sections in units of 10**(-16) cm**2
c
               dsgma(it) = pi * 2.8003d-1 * (2*fjtot+1) / xak
	       fdgn(it) =  1.d0 / ( (2*fj1js(it)+1) * (2*fj2js(it)+1) )
	    endif
	 enddo
         n12ou = nit
         do ip = 1,njsmx
            do ipp = 1,njsmx
               opc(njsmx*(ip-1)+ipp) = 0
	    enddo
	 enddo
c
c     loop over total parity
c
	 do i = 1,2
c
c     determine the sets of quantum numbers associated with the basis
c     states in both the space-fixed and body-fixed representations
c
	    call bases (fj12  ,fj12mx,pprt  ,par(i),fjtot  ,nj12
     &                 ,idxbf ,idxsf ,fl    ,xmg   ,xnorm  ,nv
     &                 ,nblk  ,ndmblk)
            if (nv .ne. 0) then
c
c     calculate all algebraic coefficients required for the evaluation of
c     the potential and for the body fixed to space fixed transformation
c     (see Launay 1977)
c
	       call algeb (fj12  ,fjtot ,par(i),idxbf ,idxsf ,xmg
     &                    ,fl    ,xnorm ,nv    ,nq12  ,q1q   ,q2q
     &                    ,q12   ,q1m   ,q2m   ,fmu   ,xred  ,pprt
     &                    ,nj12  ,glm   ,nblk  ,ndmblk,ndbmat,cj
     &                    ,cr    ,ntrt  ,ir1p  ,ir2p  ,nvou  ,idbglm
     &                    ,enbf  ,ak    ,akr   ,ekjj)
c
c     integration of coupled equations
c
	       if (integrator .le. 0) then
		  call dvglr (aa    ,bb    ,cc    ,dd    ,ee    ,nvou
     &                       ,nv    ,nq12  ,nnq12 ,nptpot,rpot  ,vpot
     &                       ,cm    ,vv    ,glm   ,ndbmat,idbglm,ndmblk
     &                       ,cr    ,enbf  ,nblk  ,ntrt  ,ir1p  ,ir2p
     &                       ,ndbpot,idxbf ,vmat  ,ipiv  ,work  ,rmatch)
                  call mtcopy (zz,bb,nv,nv,nv,nv)
               else
                  initjo = 1
		  call jonsom (zz    ,aa    ,bb    ,cc    ,dd    ,ee
     &                        ,nv    ,rdb   ,rfn   ,rmatch,fpt   ,nv
     &                        ,initjo,itest ,integrator,nblk
     &                        ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,glm
     &                        ,vv    ,nq12  ,nnq12 ,nptpot,rpot  ,vpot
     &                        ,cm    ,ndbpot,ir1p  ,ir2p  ,idxbf ,ntrt)
               endif
c
	       call bf2sf (aa,zz,nv,idxbf,idxsf,cj,nv  )
c
c     calculation of K and S matrices
c     note that the K (reactance) matrix is real and symmetric,
c     the S (scattering) matrix is unitary and symmetric; the t
c     (transmission) matrix is related to the S matrix through
c     the matrix relation T = 1 - S
c
	       call kmatrx (aa    ,zz    ,cc    ,rmatch,nv    ,ak
     &                     ,akr   ,fl    ,sb    ,sbp   ,sn    ,snp
     &                     ,ipiv  ,work)
	       call smatrx (aa    ,zz    ,cc    ,nvou
     &                     ,nv    ,2     ,it2)
c
c     prints K matrix
c
               if (ikmat .ge. 1)
     &            call mtlook (cc,nvou,1,nvou,1,nv,8,2,'k-matrix')
c
c     prints |T|**2 matrix
c
	       if (it2 .ge. 1)
     &            call mtlook (cc,nvou,1,nvou,1,nv,8,2,'t**2-mat')
c
c     contracts T matrix
c
               ips = 0
               do it = 1,n12ou
                  isom = 0
                  ndbsom(it) = ips
                  do iv = 1,nvou
                     iw = idxsf(iv)
                     if (dabs(fn1js(it)-fn1jj(iw)) .gt. eps) go to 250
                     if (dabs(fl1js(it)-fl1jj(iw)) .gt. eps) go to 250
                     if (dabs(fs1js(it)-fs1jj(iw)) .gt. eps) go to 250
                     if (dabs(fj1js(it)-fj1jj(iw)) .gt. eps) go to 250
                     if (dabs(fj2js(it)-fj2jj(iw)) .gt. eps) go to 250
                     if (dabs(fs2js(it)-fs2jj(iw)) .gt. eps) go to 250
                     if (dabs(fl2js(it)-fl2jj(iw)) .gt. eps) go to 250
                     if (dabs(fn2js(it)-fn2jj(iw)) .gt. eps) go to 250
                     isom = isom + 1
                     ips = ips + 1
                     ivv(ips) = iv
 250                 continue
		  enddo
                  nsom(it) = isom
	       enddo
c
c     check if dimension exceeded
c
               if (n12ou .gt. nvmx) write (6, 9420) n12ou,nvmx
               do ip = 1,nvou
                  do it = 1,n12ou
		     sum = 0.d0
                     n = nsom(it)
                     if (n .ne. 0) then
                        ndb = ndbsom(it)
                        do iq = 1,n
                           iv = ivv(ndb+iq)
			   sum = sum + cc(nv*(iv-1)+ip)
			enddo
		     endif
		     ee(nvmx*(it-1)+ip) = sum
		  enddo
	       enddo
               do ip = 1,n12ou
                  do it = 1,n12ou
		     sum = 0.d0
                     n = nsom(it)
                     if (n .ne. 0) then
                        ndb = ndbsom(it)
                        do iq = 1,n
                           iv = ivv(ndb+iq)
			   sum = sum + ee(nvmx*(ip-1)+iv)
			enddo
		     endif
		     dd (nvmx*(ip-1)+it) = sum
		  enddo
	       enddo
c
c     computes partial probabilities
c
               do ip = 1,n12ou
                  do ipp = 1,n12ou
                     opc(njsmx*(ip-1)+ipp) = opc(njsmx*(ip-1)+ipp)
     &                                     + dd(nvmx*(ip-1)+ipp)
		  enddo
	       enddo
	    endif
	 enddo
         do ip = 1,n12ou
            do ipp = 1,n12ou
	       opc(njsmx*(ip-1)+ipp) = opc(njsmx*(ip-1)+ipp)
     &                               * fdgn(ip)
	    enddo
	 enddo
         opcmx = 0
         do ip = 1,n12ou
            do ipp = 1,n12ou
               if (ip .ne. ipp) then
                  opcc = opc(njsmx*(ip-1)+ipp)
                  if (opcc .gt. opcmx)  opcmx = opcc
	       endif
	    enddo
	 enddo
c
         if (iopc .ge. 1) then
	    call mtlook (opc,n12ou,1,n12ou,1,njsmx,8,2,'opacity ')
         endif
c
c     computes partial cross sections
c
         do ip = 1,n12ou
            do ipp = 1,n12ou
               opc(njsmx*(ip-1)+ipp) = opc(njsmx*(ip-1)+ipp)
     &                               * dsgma(ip)
	    enddo
	 enddo
c
         if (isctp .ge. 1) then
	    call mtlook (opc,n12ou,1,n12ou,1,njsmx,8,2,'part-xs ')
         endif
c
c     calcul des sections efficaces totales
c
         do ip = 1,n12ou
            do ipp = 1,n12ou
	       sct(njsmx*(ip-1)+ipp) = sct(njsmx*(ip-1)+ipp)
     &                               + opc(njsmx*(ip-1)+ipp)
     &                               * fjps
	    enddo
	 enddo
c
         if (isctot .ge. 1) then
	    call mtlook (sct,n12ou,1,n12ou,1,njsmx,8,2,'tot-xs  ')
         endif
c
c   computes total inelasticities and energy transfers
c   (skipped if itr = 0)
c
         if (itr .ge. 1) then
            do ipp = 1,n12ou
               i1p = idx1rg(ipp)
               i2p = idx2rg(ipp)
               xa = 0
               xb = 0
               xc = 0
               xd = 0
               do ip = 1,n12ou
                  i1 = idx1rg(ip)
                  i2 = idx2rg(ip)
                  dele = (e1k(i1) + e2k(i2) - e1k(i1p) - e2k(i2p))
     &                 / dmskua
                  xa = xa + opc(njsmx*(ipp-1)+ip) * dele
                  xb = xb + sct(njsmx*(ipp-1)+ip) * dele
                  xc = xc + opc(njsmx*(ipp-1)+ip)
                  xd = xd + sct(njsmx*(ipp-1)+ip)
	       enddo
	       sctinp(ipp)  = xc - opc(njsmx*(ipp-1)+ipp)
               sctint (ipp) = xd - sct(njsmx*(ipp-1)+ipp)
               trp (ipp) = xa / sctinp(ipp)
	       trt(ipp)  = xb / sctint(ipp)
	    enddo
            write (6, 9510)
c
c     print the partial and total inelasticities (10**(-16) cm**2)
c
            write (6, 9500) (sctinp(ipp),ipp=1,n12ou)
            write (6, 9500) (sctint(ipp),ipp=1,n12ou)
            write (6, 9520)
c
c     print the partial and total energy transfers (k)
c
            write (6, 9500) (trp   (ipp),ipp=1,n12ou)
            write (6, 9500) (trt   (ipp),ipp=1,n12ou)
	 endif
c
c     compare current and maximum values of the partial probabilities
c     (test of convergence of the partial wave expansion)
c     and of the total angular momentum
c
         if (opcmx .ge. opcmn) then
            if (fjtot .lt. (fjfn-eps)) go to 200
	 endif
      call mtlook (sct,n12ou,1,n12ou,1,njsmx,8,2,'tot-xs  ')
      enddo
      return
c
c    error message when dimension exceeded
c
 9000 format (' *** WARNING NNQ12 = ',i5,' GREATER THAN NQ12MX = ',i5)
 9010 format (' **** WARNING NPT = ',i5,' GREATER THAN NPTMX = ',i5)
 9020 format (' **** WARNING NQ12MX = ',i5
     &       ,' GREATER THAN  150: REDIMENSION COMMON/VIB/')
 9030 format (' **** WARNING NJ12MX = ',i5
     &       ,' GREATER THAN  90: REDIMENSION COMMON/VIB/')
 9050 format (' **** WARNING N1MX = ',i5
     &       ,' GREATER THAN  30: REDIMENSION E1K ETC. IN PRNCP')
 9060 format (' **** WARNING N2MX = ',i5
     &       ,' GREATER THAN  30: REDIMENSION E2K ETC. IN PRNCP')
 9070 format (' **** WARNING NVMX = ',i5
     &       ,' GREATER THAN  801: REDIMENSION TRP ETC. IN PRNCP')
 9200 format (i3/(6d13.6))
 9210 format (i3/(10d13.6))
 9300 format (10f5.2)
 9310 format (' ',13f10.3)
 9320 format (5x,'   N1        L1        S1        J1     '
     &    ,'   N2        L2        S2        J2        J12    ENERGY ')
 9350 format (5x,'   N1        N2        Q1        Q2        Q12'
     &       ,'       Q1        Q2        MU ')
 9400 format (' ***** ENERGY  = ',f15.4,' DEGREES KELVINS **********')
 9410 format (30x,'********** PARTIAL WAVE J = ',f5.1,' ********')
 9420 format (' **** WARNING N12OU = ',i5,' GREATER THAN NVMX  = ',i5)
 9500 format (' ',1p12d10.3)
 9510 format ('   PARTIAL AND TOTAL INELASTICITIES  ')
 9520 format ('   PARTIAL AND TOTAL ENERGY TRANSFERS ')
      end
c-----------------------------------------------------------------------
      subroutine dvglr (aa    ,bb    ,cc    ,dd    ,ee    ,nvou
     &                 ,nv    ,nq12  ,nnq12 ,nptpot,rpot  ,vpot
     &                 ,cm    ,vv    ,glm   ,ndbmat,idbglm,ndmblk
     &                 ,cr    ,enbf  ,nblk  ,ntrt  ,ir1p  ,ir2p
     &                 ,ndbpot,idxbf ,vmat  ,ipiv  ,work  ,rend)
c
c     dvglr is called from prncp
c     integration of the second order differential equations using
c     the method of De Vogelaere (see Launay 1976)
c     integration starts in the classically forbidden region, at
c     intermolecular distance r = rdb and continues to rfn, where
c     the asymptotic form of the solutions should have been attained
c
c     calls cpot to calculate the potential energy curves at the
c     current value of r
c
c     calls vblk to calculate the coupling matrix elements,
c     incorporating the interaction potential, centrifugal potential,
c     and relative collision energy
c
      implicit double precision (a-h,o-z)
c
c     integration parameters
c
      common /intgr/ rdb,rfn,fpt,enfdml,rmin,stbmx,fllmx,enmin,nstbmx
      common /steps/ step,eneff,enstb,deltae
c
c     printout parameters
c
      common /ecrit/ ikmat,it2,iopc,isctp,isctot,itest,iradin
      dimension
     &     aa    (*),bb    (*),cc   (*),dd    (*),ee    (*),nptpot(*)
     &    ,rpot  (*),vpot  (*),cm   (*),vv    (*),glm   (*),ndbmat(*)
     &    ,idbglm(*),ndmblk(*),cr   (*),enbf  (*),ir1p  (*),ir2p  (*)
     &    ,ndbpot(*),idxbf (*),vmat  (*)
      dimension ipiv(1:nv),work(1:nv)
c
      npas = 0
      scstb = 0.d0
      nstblz = 0
      ncstb = 0
      ncps = 0
      del = 1.d0
      call stepper (rdb,hp)
      h = hp
      hh = hp / 2.d0
      rr = rdb
      delhh = -hh * del
      r = rr + hh
      do i = 1,nv*nv
	 aa(i) = 0.
	 bb(i) = 0.
	 cc(i) = 0.
	 dd(i) = 0.
	 ee(i) = 0.
      enddo
      call cpot (nnq12 ,nptpot,rpot  ,vpot  ,cm    ,rr
     &          ,vv    ,ndbpot,eneff ,enstb)
      do ib = 1,nblk
	 call vblk (ib    ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,idxbf
     &             ,glm   ,vmat  ,vv    ,rr    ,nq12  ,nnq12 ,nv)
         nn = ndbmat(ib)
         ndm = ndmblk(ib)
         ki = nn
         imat = 1
         do i = 1,ndm
            ki = ki + 1
            dd(ki+(ki-1)*nv) = vmat(imat) * delhh
            imat = imat + ndm + 1
	 enddo
      enddo
      do i = 1,nv
         bb(i+(i-1)*nv) = del
      enddo
  100 continue
      npas = npas + 1
      ncps = ncps + 1
      hh = hp / 2.
      hpcar = hp * hp
      h1 = hpcar * (3.d0 + hp/h) / 24.d0
      h2 = - hpcar * hp / 24.d0 / h
      h3= hpcar / 3.d0
      h4 = hp / 6.d0
      h5 = 4.d0 * h4
      h6 = hpcar / 6.d0
      do i = 1,nv*nv
         cc(i) = aa(i) + hh * bb(i) + h1 * ee(i) + h2 * dd(i)
      enddo
      rr = r + hh
      call cpot (nnq12 ,nptpot,rpot  ,vpot  ,cm    ,rr
     &          ,vv    ,ndbpot,eneff ,enstb)
      do ib = 1,nblk
	 call vblk (ib    ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,idxbf
     &             ,glm   ,vmat  ,vv    ,rr    ,nq12  ,nnq12 ,nv)
         nn = ndbmat(ib)
         ndm = ndmblk(ib)
	 call fmmult (dd(nn+1),vmat,cc(nn+1)
     &               ,ndm,ndm,nv ,nv,ndm,nv)
      enddo
      if (ntrt .ne. 0) then
         rr2 = rr * rr
         do i = 1,ntrt
            i1p = ir1p(i)
            i2p = ir2p(i)
            crm = cr(i) / rr2
            crp = crm
	    do j = 0,nv*(nv-1),nv
               dd(i2p+j) = dd(i2p+j) + crm * cc(i1p+j)
               dd(i1p+j) = dd(i1p+j) + crp * cc(i2p+j)
	    enddo
	 enddo
      endif
      do i = 1,nv*nv
         aa(i) = aa(i) + hp * bb(i) + h6 * ee(i) + h3 * dd(i)
      enddo
c
      r = rr + hh
      call cpot (nnq12 ,nptpot,rpot  ,vpot  ,cm    ,r
     &          ,vv    ,ndbpot,eneff ,enstb)
      do ib = 1,nblk
         nn = ndbmat(ib)
	 call vblk (ib    ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,idxbf
     &             ,glm   ,vmat  ,vv    ,r     ,nq12  ,nnq12 ,nv)
         ndm = ndmblk(ib)
	 call fmmult (cc(nn+1),vmat,aa(nn+1)
     &               ,ndm,ndm,nv ,nv,ndm,nv)
      enddo
      if (ntrt .ne. 0) then
         r2 = r * r
         do i = 1,ntrt
            i1p = ir1p(i)
            i2p = ir2p(i)
            crm = cr(i) / r2
            crp = crm
	    do j = 0,nv*(nv-1),nv
               cc(i2p+j) = cc(i2p+j) + crm * aa(i1p+j)
               cc(i1p+j) = cc(i1p+j) + crp * aa(i2p+j)
	    enddo
	 enddo
      endif
c
c     computes integration step
c
      h = hp
      call stepper (r,hp)
      do i = 1,nv*nv
         bb(i) = bb(i) + h4 * ee(i) + h5 * dd(i) + h4 * cc(i)
         ch = cc(i)
         cc(i) = ee(i)
         ee(i) = ch
      enddo
c
c     stabilisation
c
      akstb = 0.
      if (enstb .lt. 0.) akstb = dsqrt(-enstb)
      scstb = scstb + hp * akstb
      ncstb = ncstb + 1
      if (ncstb .ge. nstbmx .or. scstb .ge. stbmx) then
	 call stblz (aa,bb,cc,dd,ee,nv,ipiv,work)
         nstblz = nstblz + 1
         scstb = 0.d0
         ncstb = 0
      endif
      if (r .lt. rfn) go to 100
      rend = r
      call stblz (aa,bb,cc,dd,ee,nv,ipiv,work)
      if (itest .ge. 1) then
c
c     optional print out: collision energies (relative to individual
c     channels); number of radial integration points; number of
c     stabilizations
c
c     write (6, 6001) (enbf(i),i=1,nv)
c 6001 format(' ',1p12d10.3)
c
         write (6, 9000) npas, nstblz
      endif
      return
 9000 format (' ',i4,' INTEGRATION POINTS ',i4,' STABILISATIONS')
      end
c-----------------------------------------------------------------------
      subroutine vblk (iblk  ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,idxbf
     &                ,glm   ,vmat  ,vv    ,r     ,nq12  ,nnq12 ,nv)
c
c     vblk is called from dvglr
c     calculates the potential matrix elements, incorporating the
c     contributions of the interaction potential, the centrifugal
c     potential, and the relative collision energy
c
      implicit double precision (a-h,o-z)
      common /vib/ fn2q(150),fnp2q(150),fn1jj(90),fn2jj(90),fj2jj(90)
      dimension idbglm(*),ndmblk(*),ndbmat(*),cr(*),enbf(*),glm(*)
     &         ,vmat  (*),vv    (*),idxbf (*)
c
      np = idbglm(iblk)
      ndm = ndmblk(iblk)
      do i = 1,ndm
         do ip = 1,i
	    sum = 0.d0
            icon = 0
            do iq = 1,nnq12
               if (icon .lt. nq12) then
c
c     check the supplementary quantum number of system 2
c     the corresponding potential curve must then be used
c
		  kii = idxbf(ndbmat(iblk)+i)
		  kip = idxbf(ndbmat(iblk)+ip)
		  if (((dabs(fn2jj(kii)-fn2q (iq)) .le. 0.01d0) .and.
     &                 (dabs(fn2jj(kip)-fnp2q(iq)) .le. 0.01d0)) .or.
     &                ((dabs(fn2jj(kip)-fn2q (iq)) .le. 0.01d0) .and.
     &                 (dabs(fn2jj(kii)-fnp2q(iq)) .le. 0.01d0))) then
c
                     np = np + 1
c
c      combine algebraic coefficients glm with the appropriate
c      intermolecular-distance-dependent potential energy curve v
c
		     sum = sum + glm(np) * vv(iq)
                     icon = icon + 1
		  endif
	       endif
	    enddo
	    vmat(ndm*(i-1)+ip) = sum
	    vmat(ndm*(ip-1)+i) = sum
	 enddo
      enddo
c
      ien = ndbmat(iblk)
      icr = ndbmat(iblk) + nv
      r2 = r * r
      imat = 1
      do i = 1,ndm
         ien = ien + 1
         icr = icr + 1
c
c     combine the contributions of the interaction potential, the
c     centrifugal potential, and the relative barycentric energy
c
         vmat(imat) = vmat(imat) + cr(icr) / r2 - enbf(ien)
         imat = imat + ndm + 1
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine algeb (fj12  ,fjtot ,par   ,idxbf ,idxsf ,xmg
     &                 ,fl    ,xnorm ,nv    ,nq12  ,q1q   ,q2q
     &                 ,q12   ,q1m   ,q2m   ,fmu   ,xred  ,pprt
     &                 ,nj12  ,glm   ,nblk  ,ndmblk,ndbmat,cj
     &                 ,cr    ,ntrt  ,ir1p  ,ir2p  ,nvou  ,idbglm
     &                 ,enbf  ,ak    ,akr   ,ekjj)
c
c     algeb is called from dvglr
c     calculates all the algebraic coefficients necessary to evaluate the
c     potential matrix elements and to perform the body-fixed to space-
c     fixed transformation (Launay 1976, 1977)
c
      implicit double precision (a-h,o-z)
      common /intgr/ rdb,rfn,fpt,enfdml,rmin,stbmx,fllmx,enmin,nstbmx
      common /ecrit/ ikmat,it2,iopc,isctp,isctot,itest,iradin
      common /ndim/ nptmx,nvmx,njsmx,nj12mx,nq12mx,ndgnmx,ncjmx,nglmmx
     &             ,nblkmx,n1mx,n2mx
      common /vib/ fn2q(150),fnp2q(150),fn1jj(90),fn2jj(90),fj2jj(90)
      dimension
     &     fj12  (*),idxbf (*),idxsf (*),xmg   (*),fl    (*),xnorm (*)
     &    ,q1q   (*),q2q   (*),q12   (*),q1m   (*),q2m   (*),fmu   (*)
     &    ,xred  (*),pprt  (*),glm   (*),ndmblk(*),ndbmat(*),cj    (*)
     &    ,cr    (*),ir1p  (*),ir2p  (*),idbglm(*),enbf  (*),ak    (*)
     &    ,akr   (*),ekjj  (*)
      data eps /0.01d0/
      xl(xxx) = 2.*xxx+1.
      root2 = dsqrt(2.d0)
c
c     computes coefficients of BF --> SF unitary transformation
c
      icj = 0
      do isf = 1,nv
         ijjsf = idxsf(isf)
         do ibf = 1,nv
            ijjbf = idxbf(ibf)
            if (ijjbf .eq. ijjsf) then
               nx = fjtot - xmg(ibf) + eps
	       sgn = 1.d0
               if (mod(nx,2) .gt. 0) sgn = -sgn
               icj = icj + 1
	       cj(icj) = sgn
     &                 * dsqrt(xl(fl(isf)))
     &                 * f3j (fj12(ijjbf),fjtot,fl(isf)
     &                       ,xmg(ibf),-xmg(ibf),0.d0)
     &                 / xnorm(ibf)
     &                 * root2
	    endif
	 enddo
      enddo
      ncj = icj
c
c     check if input dimension exceeded
c
      if (ncj .gt. ncjmx) write (6, 9000) ncj,ncjmx
c
c     computes algebraic coefficients for interaction potential
c     (incorporate the reduced matrix elements xred evaluated in prncp)
c
      iglm = 0
      do i = 1,nv
         do ip = 1,i
            if (dabs(xmg(i)-xmg(ip)) .le. eps) then
               ii = idxbf(i)
               iip = idxbf(ip)
               ixp = (ii * (ii-1) / 2 + iip - 1) * nq12
               do iqm = 1,nq12
                  iglm = iglm + 1
		  sum = 0.d0
                  nx = q1m(iqm) + q2m(iqm) + fj12(ii) - xmg(i) + eps
		  sgn = 1.d0
                  if (mod(nx,2) .gt. 0) sgn = -sgn
		  aaa = root2
                  fmum = fmu(iqm)
                  q1mm = q1m(iqm)
		  if (dabs(fmum) .lt. eps) aaa = 1.d0
                  q2mm = q2m(iqm)
                  do iqq = 1,nq12
                     q1qq = q1q(iqq)
                     q2qq = q2q(iqq)
                     if (dabs(q1mm-q1qq) .le. eps .and.
     &                   dabs(q2mm-q2qq) .le. eps)
     &                   sum = sum
     &                       + sgn
     &                       * aaa
     &                       * f3j (q1qq,q2qq,q12(iqq)
     &                             ,-fmum,fmum,0.d0)
     &                       * f3j (fj12(ii),q12(iqq),fj12(iip)
     &                             ,-xmg(i),0.d0,xmg(i))
     &                       * xred (ixp+iqq)
		  enddo
c
c     store the algebraic coefficients for subsequent use
c
		  glm(iglm) = sum
	       enddo
	    endif
	 enddo
      enddo
      nglm = iglm
c
c     check if input dimension exceeded
c
      if (nglm .gt. nglmmx) write (6, 9010) nglm,nglmmx
      ndbmat(1) = 0
      do i = 1,nblk
	 ndbmat(i+1) = ndbmat(i) + ndmblk(i)
      enddo
      np = 0
      do ib = 1,nblk
         ndm = ndmblk(ib)
         idbglm(ib) = np
         do i = 1,ndm
	    np = np + i * nq12
	 enddo
      enddo
c
c     calculation of l**2 matrix elements
c
      icm = 0
      ic = 0
      do i1 = 1,nv
         do i2 = i1 , nv
	    sum = 0.d0
            ii1 = idxbf(i1)
            ii2 = idxbf(i2)
            f12 = fj12(ii1)
            if (ii1 .eq. ii2) then
               xmg1 = xmg(i1)
               xmg2 = xmg(i2)
               if (dabs(xmg1-xmg2) .le. eps) then
                  ic = ic + 1
		  sum  = fjtot * (fjtot+1)
     &                 + f12 * (f12+1)
     &                 - 2 * xmg1 * xmg2
                  if (dabs(xmg1-0.5d0) .le. eps) then
                     nx = fjtot + 3 * f12 + eps
		     sgn = 1.d0
                     if (mod(nx,2) .gt. 0) sgn = -sgn
		     sum = sum
     &                   - sgn
     &                   * par * pprt(ii1)
     &                   * dsqrt ( (f12*(f12+1)     + 0.25d0)
     &                           * (fjtot*(fjtot+1) + 0.25d0) )
		  endif
		  cr(ic+nv) = sum
	       endif
               if (dabs(xmg2-xmg1-1) .le. eps) then
                  icm = icm + 1
		  sum = - dsqrt( ( fjtot*(fjtot+1) - xmg1*xmg2)
     &                           * (f12*(f12+1)    - xmg1*xmg2) )
		  if (dabs(xmg1) .lt. eps) sum = sum * root2
		  cr(icm) = sum
                  ir1p(icm) = i1
                  ir2p(icm) = i2
	       endif
	    endif
	 enddo
      enddo
      ntrt = icm
c
c     definition of energy and wave vectors for each channel
c
      nvou = 0
      enmin = 1.d72
      do i = 1,nv
         ibf = idxbf(i)
         isf = idxsf(i)
c
c     relative collision energy
c
         enbf(i) = enfdml - (ekjj(ibf)-ekjj(1))
         ww = enfdml - (ekjj(isf)-ekjj(1))
         if (enmin .gt. ww) enmin = ww
         ak(i) = dsqrt(dabs(ww))
         akr(i) = dsqrt(ak(i))
c
c     count the number of open channels
c
         if (ww .gt. 0.0d0) nvou = nvou+1
         if (ww .lt. 0.0d0) ak(i) = -ak(i)
      enddo
c
c     optional print out of the numbers of channels, body-fixed to
c     space-fixed transformation coefficients, and potential coefficients
c
      if (itest .gt. 0) write (6, 9020) nv,ncj,nglm
c
      return
 9000 format (' ***** WARNING NCJ = ',i5,'GREATER THAN NCJMX = ',i5)
 9010 format (' ***** WARNING NGLM = ',i5,'GREATER THAN NGLMMX = ',i5)
 9020 format (' NV = ',i7,5x,'NCJ = ',i7,5x,'NGLM = ',i7)
      end
c-----------------------------------------------------------------------
      subroutine bases (fj12  ,fj12mx,pprt  ,par   ,fjtot ,nj12
     &                 ,idxbf ,idxsf ,fl    ,xmg   ,xnorm ,nv
     &                 ,nblk  ,ndmblk)
c
c     this routine computes quantum numbers
c     omega(xmg) et l(fl)  of SF and BF representations
c     computes also channel number and other quantum numbers
c       using idxbf and idxsf arrays
c
      implicit double precision (a-h,o-z)
      common /intgr/ rdb,rfn,fpt,enfdml,rmin,stbmx,fllmx,enmin,nstbmx
      common /ndim/  nptmx,nvmx,njsmx,nj12mx,nq12mx,ndgnmx,ncjmx,nglmmx
     &              ,nblkmx,n1mx,n2mx
      dimension fj12  (*),pprt  (*),idxbf (*),idxsf (*),fl    (*)
     &         ,xmg (*)  ,xnorm (*),ndmblk(*)
      data eps /0.01d0 /
      root2 = dsqrt(2.d0)
      i = fjtot + eps
      xmgmn = 0.
      if (dabs(fjtot-i) .gt. eps) xmgmn = 0.5d0
      xmgmx = dmin1(fjtot,fj12mx)
      nv = 0
      nblk = 0
      xox = xmgmn
 20   continue
      ndm = 0
      do i = 1,nj12
	 if (fj12(i) .ge. (xox  - eps)) then
            nx = fjtot + 3 * fj12(i) + eps
            pr1 = par
            if (mod(nx,2) .gt. 0) pr1 = -par
	    if ((dabs(xox) .ge. eps) .or. (pr1*pprt(i) .ge. 0.d0)) then
               ndm = ndm + 1
               nv = nv + 1
	       xnorm(nv) = 1.d0
	       if ((dabs(xox) .lt. eps) .and. (pr1*pprt(i) .gt. 0.d0))
     &              xnorm(nv) = root2
               idxbf(nv) = i
	       xmg(nv) = xox
	    endif
	 endif
      enddo
      if (ndm .ne. 0) then
         nblk = nblk + 1
         ndmblk(nblk) = ndm
      endif
      xox =  xox + 1
      if (xox .lt. (xmgmx+eps)) go to 20
c
c     check if dimension exceeded ; nblkmx is defined in dime
c
      if (nblk .gt. nblkmx) write (6, 9000) nblk,nblkmx
      if (nv .eq. 0) return
      nv = 0
      fllmx = 0.d0
      do i = 1,nj12
         flmn = dabs(fj12(i)-fjtot)
         flmx = fj12(i) + fjtot
         flmxu = flmx * (flmx+1.)
         if (fllmx .lt. flmxu) fllmx = flmxu
	 flx = flmn
 45      continue
	 nx = flx + eps
         pr1 = par
         if (mod(nx,2) .gt. 0) pr1 = -par
         if (pr1*pprt(i) .ge. 0.d0) then
            nv = nv + 1
	    fl(nv) = flx
            idxsf(nv) = i
	 endif
	 flx = flx + 1
	 if (flx .lt. (flmx+eps)) go to 45
      enddo
c
c     check if dimension exceeded
c
      if (nv .gt. nvmx) write (6, 9010) nv,nvmx
 9000 format (' **** WARNING NBLK = ',i5,' GREATER THAN NBLKMX = ',i5)
 9010 format (' ***** WARNING NV = ',i5,' GREATER THAN NVMX = ',i5)
      return
      end
c-----------------------------------------------------------------------
      subroutine kmatrx (aa    ,bb    ,rmat  ,r     ,nv    ,ak
     &                  ,akr   ,fl    ,sb    ,sbp   ,sn    ,snp
     &                  ,ipiv  ,work)
c
c     calculates the k matrix elements
c
      implicit double precision (a-h,o-z)
      dimension aa(nv,nv),bb(nv,nv),rmat(nv,nv),ak(*),akr(*),fl(*)
     &         ,sb(*),sbp(*),sn(*),snp(*),ipiv(1:nv),work(1:nv)
c
      do i = 1,nv
	 aki = ak(i)
         rcqm = akr(i)
	 rcqm = 1.d0 / rcqm
	 xx = aki * r
	 aki1 = aki
	 if (aki1 .lt. 0.) aki1 = -aki
	 aai = fl(i)*(fl(i)+1.d0)
	 call besph2 (bj,bjp,bn,bnp,aai,xx,1,0)
	 sb(i)  = rcqm * bj
	 sbp(i) = rcqm * aki1 * bjp
	 sn(i)  = rcqm * bn
	 snp(i) = rcqm * aki1 * bnp
      enddo
      do j = 1,nv
         do i  = 1,nv
            aa(i,j) = bb(i,j) * sb(j)
            bb(i,j) = bb(i,j) * sn(j)
	 enddo
      enddo
      do i = 1,nv
         bb(i,i) = bb(i,i) - snp(i)
         aa(i,i) = aa(i,i) - sbp(i)
      enddo
c
      call dgetrf (nv,nv,bb,nv,ipiv,info)
      call dgetri (nv,bb,nv,ipiv,work,nv,info)
c
      call fmmult (rmat,bb,aa, nv,nv,nv, nv,nv,nv)
      return
      end
c-----------------------------------------------------------------------
      subroutine stblz (aa,bb,cc,dd,ee,nv,ipiv,work)
c
c     stblz is called from the de vogelaere integration routine dvglr
c     performs the stabilization of the solutions to
c     the coupled second order differential equations
c
      implicit double precision (a-h,o-z)
      dimension aa(nv,*),bb(nv,*),cc(nv,*),dd(nv,*),ee(nv,*)
      dimension ipiv(1:nv),work(1:nv)
c
      call dgetrf (nv,nv,aa,nv,ipiv,info)
      call dgetri (nv,aa,nv,ipiv,work,nv,info)
c
      call fmmult (cc,bb,aa, nv,nv,nv, nv,nv,nv)
      call mtcopy (bb,cc, nv,nv, nv,nv)
c
      call fmmult (cc,dd,aa, nv,nv,nv, nv,nv,nv)
      call mtcopy (dd,cc, nv,nv, nv,nv)
c
      call fmmult (cc,ee,aa, nv,nv,nv, nv,nv,nv)
      call mtcopy (ee,cc, nv,nv, nv,nv)
c
      do i1 = 1,nv
         do i2 = 1,nv
            aa(i1,i2) = 0.d0
	 enddo
      enddo
      do i = 1,nv
         aa(i,i) = 1.d0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine smatrx (ure   ,uim   ,usq   ,nop   ,nmx   ,keyout,ibug)
c#######################################################################
c#    input  : k-matrix or eigenvectors and eigenphases of k-matrix.   #
c#    output : s-matrix or t-matrix.                                   #
c#---------------------------------------------------------------------#
c#    s = (1+i*k)/(1-i*k) , t = 1-s                                    #
c#    re(s) = (1-k*k)/(1+k*k)                 im(s) = 2*i*k/(1+k*k)    #
c#---------------------------------------------------------------------#
c#    ure    : real part of the S or T matrix                          #
c#    uim    : imaginary part of the S or T matrix                     #
c#    usq    : on input,  contains the K-matrix                        #
c#             on output, contains the norms squared of S or T matrix  #
c#                        elements                                     #
c#    nop    : number of channels                                      #
c#    nmx    : dimension of arrays in calling program                  #
c#    keyout : le 1 for s-matrix, otherwise t-matrix                   #
c#    ibug   : positive for printing the matrix                        #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (ncmx1=801,lwork=64*ncmx1)
      dimension ure(nmx,*),uim(nmx,*),usq(nmx,*)
      dimension ipiv(ncmx1),work(lwork)
c
      call fmmult (ure,usq,usq, nop,nop,nop ,nmx,nmx,nmx)
      do i = 1,nop
	 ure(i,i) = 1.d0 + ure(i,i)
      enddo
      do i = 1,nop
	 do j = 1,nop
	    ure(i,j) = 0.5d0 * ure(i,j)
	 enddo
      enddo
c
      call dgetrf (nop,nop,ure,nmx,ipiv,info)
      call dgetri (nop,ure,nmx,ipiv,work,lwork,info)
c
      call fmmult (uim,ure,usq,nop,nop,nop,nmx,nmx,nmx)
      do i = 1,nop
	 ure(i,i) = -1.d0 + ure(i,i)
      enddo
c
      if (keyout .ge. 2) then
c
	 do i = 1,nop
	    do j = 1,nop
	       ure(i,j) = -ure(i,j)
	       uim(i,j) = -uim(i,j)
	    enddo
	 enddo
	 do i = 1,nop
	    ure(i,i) = 1.d0 + ure(i,i)
	 enddo
c
      endif
c
      do i = 1,nop
	 do j = 1,nop
	    usq(i,j) = ure(i,j)**2 + uim(i,j)**2
	 enddo
      enddo
c
      if (ibug .ge. 1) then
	 if (keyout .eq. 1) call matprt (usq,nop,nop,nmx,8,1,'|s|**2  ')
	 if (keyout .ne. 1) call matprt (usq,nop,nop,nmx,8,2,'|t|**2  ')
      endif
      if (ibug .ge. 2) then
	 if (keyout .eq. 1) call matprt (ure,nop,nop,nmx,8,1,'re(s)   ')
	 if (keyout .eq. 1) call matprt (uim,nop,nop,nmx,8,1,'im(s)   ')
	 if (keyout .ne. 1) call matprt (ure,nop,nop,nmx,8,1,'re(t)   ')
	 if (keyout .ne. 1) call matprt (uim,nop,nop,nmx,8,1,'im(t)   ')
      endif
      return
      end
c -------------------------------------------------------------------------
      subroutine bf2sf (aa,bb,nv,idxbf,idxsf,cj,ncmx)
c
c     transformation of the solution from the BF to the SF frame
c
      implicit double precision (a-h,o-z)
      dimension aa(ncmx,*),bb(ncmx,*),idxbf(*),idxsf(*),cj(*)
      do ibf = 1,nv
	 icj = 0
	 do isf = 1,nv
	    sum = 0.
	    ijj = idxsf(isf)
	    do jbf = 1,nv
	       jjj = idxbf(jbf)
	       if (ijj .eq. jjj) then
		  icj = icj+1
		  sum = sum + bb(ibf,jbf) * cj(icj)
	       endif
	    enddo
	    aa(ibf,isf) = sum
	 enddo
      enddo
c
      do isf = 1,nv
	 icj = 0
	 do jsf = 1,nv
	    sum = 0.
	    jjj = idxsf(jsf)
	    do ibf = 1,nv
	       ijj = idxbf(ibf)
	       if (jjj .eq. ijj) then
		  icj = icj+1
		  sum = sum + cj(icj) * aa(ibf,isf)
	       endif
	    enddo
	    bb(jsf,isf) = sum
	 enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c   cubic spline interpolation routines
      subroutine spline (n,x,y,cm)
      implicit double precision ( a-h , o-z)
      common /tab/ alpha(100),beta(100),gamma(100),b(100)
      dimension x(*),y(*),cm(*)
      cm(1) = 0.d0
      cm(n) = 0.d0
      do i = 3,n
         a = x(i-1) - x(i-2)
         c = x(i) - x(i-1)
         alpha(i-2) = (a+c) / 3.d0
         beta(i-2) = c / 6.d0
         gamma(i-2) = beta(i-2)
         b(i-2) = (y(i)-y(i-1)) / c - (y(i-1)-y(i-2)) / a
      enddo
      call tridia (alpha,beta,gamma,b,cm(2),n-2)
      return
      end
c-----------------------------------------------------------------------
      subroutine tridia (alpha,beta,gamma,b,x,n)
      implicit double precision ( a-h , o-z)
      dimension alpha(*),beta(*),gamma(*),b(*),x(*)
      do i = 2,n
         rap = beta(i-1) / alpha(i-1)
         alpha(i) = alpha(i) - rap * gamma(i-1)
         b(i) = b(i) - rap * b(i-1)
      enddo
      x(n) = b(n) / alpha(n)
      do j = 2,n
         i = n - j + 1
         x(i) = (b(i) - gamma(i) * x(i+1)) / alpha(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      function spl (n,x,y,m,t)
      implicit double precision ( a-h , o-z)
      dimension x(*),y(*)
      double precision m(*)
      if (t .le. x(1)) then
         e = x(2) - x(1)
         spl = ((y(2)-y(1)) / e - m(2) * e / 6.d0) * (t - x(1)) + y(1)
      else if (t .ge. x(n)) then
         e = x(n) - x(n-1)
         spl = ((y(n)-y(n-1)) / e + m(n-1) * e / 6.d0) * (t - x(n))
     &       + y(n)
      else
         k = 1
 10      k = k + 1
         if (t .gt. x(k)) go to 10
         e = x(k) - x(k-1)
         f = x(k) - t
         g = t - x(k-1)
         spl = (m(k-1)*f*f*f + m(k)*g*g*g + (6.d0*y(k) - m(k)*e*e) * g
     &          + (6.d0*y(k-1) - m(k-1)*e*e) * f) / (6.d0*e)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine jonsom (zz    ,vv    ,yy1   ,yy4   ,qq    ,vvb
     &                  ,nch   ,rstart,rend  ,rmatch,fpt   ,ncmx
     &                  ,init  ,ibug  ,integrator,nblk
     &                  ,idbglm,ndmblk,ndbmat,cr    ,enbf  ,glm
     &                  ,vqq   ,nq12  ,nnq12 ,nptpot,rpot  ,vpot
     &                  ,cm    ,ndbpot,ir1p  ,ir2p  ,idxbf ,ntrt)
c#######################################################################
c#    integration method for second order coupled differential         #
c#    equations with variable step size.                               #
c#    b.r. johnson, j. comput. phys. (1973), 13, 445-449               #
c#    d.e. manolopoulos, j. chem. phys. (1986), 84, 6425-9             #
c#-------------------------------------------------------------------- #
c#    zz     : on input , contains log. der. matrix at rstart          #
c#           : on output, contains log. der. matrix at rend            #
c#    vv     : the matrix 2*mu*(e*i-u(x)) (dimension nch*nch)          #
c#    nch    : number of channels                                      #
c#    ncmx   : first dimension of arrays in calling program            #
c#    rstart : starting point for the integration                      #
c#    rend   : end point for the integration                           #
c#    ibug   : >  0 for information printout                           #
c#    integrator : .le. 1, johnson method                              #
c#                 .eq. 2, johnson-manolopoulos method                 #
c#                 .eq. 3, johnson method                              #
c#                          with symmetric matrix inverter             #
c#                 .ge. 4, johnson-manolopoulos method                 #
c#                          with symmetric matrix inverter             #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (ncmx1=801)
      dimension
     &     zz(ncmx,*),vv(ncmx,*),yy1(ncmx,*),yy4(ncmx,*),qq(ncmx,*)
     &    ,vvb(ncmx,*),idbglm(*),ndmblk(*),ndbmat(*),cr(*),enbf(*)
     &    ,glm(*),vqq(*),nptpot(*),rpot  (*),vpot  (*),cm   (*)
     &    ,ndbpot(*),ir1p(*),ir2p(*),idxbf(*)
      dimension y1(ncmx1),y2(ncmx1),y3(ncmx1),y4(ncmx1),vref(ncmx1)
      logical lsym
c
      rcur = rstart
      call eqnsj (vvb   ,rcur  ,nch   ,nblk  ,idbglm,ndmblk
     &           ,ndbmat,cr    ,enbf  ,glm   ,vqq   ,nq12
     &           ,nnq12 ,nptpot,rpot  ,vpot  ,cm    ,ndbpot
     &           ,ir1p  ,ir2p  ,idxbf ,ntrt  ,ncmx)
      call stepper (rcur,hcur)
      nstep = 0
      if (init .ge. 1) then
c        print *,' initialisation performed'
	 call mtinit (zz ,nch,nch,ncmx)
	 call mtinit (vv ,nch,nch,ncmx)
	 call mtinit (qq ,nch,nch,ncmx)
	 call mtinit (yy1,nch,nch,ncmx)
	 call mtinit (yy4,nch,nch,ncmx)
	 do j = 1,nch
	    zz(j,j) = dsqrt(dabs(vvb(j,j)))
c           zz(j,j) = 1.d2
	 enddo
      endif
      lsym = .false.
      if (integrator .ge. 3) lsym = .true.
c
 999  ra = rcur
      hh = hcur
      nstep = nstep + 2
c     print *,' rcur,hcur = ',rcur,hcur
      hsqo6 = hh * hh / 6.d0
      ho3 = hh / 3.d0
      foh = 4.d0 / hh
      ooh = 1.d0 / hh
c
c --- compute the reference potential and the associated propagator
c
      rc = ra + hh
      call eqnsj (vv    ,rc    ,nch   ,nblk  ,idbglm,ndmblk
     &           ,ndbmat,cr    ,enbf  ,glm   ,vqq   ,nq12
     &           ,nnq12 ,nptpot,rpot  ,vpot  ,cm    ,ndbpot
     &           ,ir1p  ,ir2p  ,idxbf ,ntrt  ,ncmx)
      if (integrator .eq. 1 .or. integrator .eq. 3) then
	 do j = 1,nch
	    vref(j) = 0.
	 enddo
      else
	 do j = 1,nch
	    vref(j) = vv(j,j)
	 enddo
      endif
      do  j = 1,nch
	 sk = dsqrt(dabs(vref(j)))
	 skh = sk * hh
	 if (vref(j) .lt. 0.) then
	    y1(j) = sk / dtan(skh)
	    y2(j) = sk / dsin(skh)
	 else if (vref(j) .gt. 0.) then
	    y1(j) = sk / dtanh(skh)
	    y2(j) = sk / dsinh(skh)
	 else
	    y1(j) = 1.d0 / hh
	    y2(j) = 1.d0 / hh
	 endif
      enddo
      do  j = 1,nch
	 y4(j) = y1(j)
	 y3(j) = y2(j)
      enddo
c
c --- calculation of y4(a,c)
c
      do j = 1,nch
	 vv(j,j) = vv(j,j)-vref(j)
      enddo
      do j = 1,nch
	 do i = 1,nch
	    vv(i,j) = -hsqo6 * vv(i,j)
	 enddo
	 vv(j,j) = 1.d0 + vv(j,j)
      enddo
      call invers (vv,qq,nch,ncmx,lsym)
      do j = 1,nch
	 do i = 1,nch
	    qq(i,j) = foh * qq(i,j)
	 enddo
	 qq(j,j) = qq(j,j) - foh
      enddo
      do j = 1,nch
	 do i = 1,nch
	    yy4(i,j) = qq(i,j)
	 enddo
	 yy4(j,j) = y4(j) + yy4(j,j)
      enddo
c
c --- calculation of y1(a,c)
c
      call mtcopy (vv,vvb,nch,nch,ncmx,ncmx)
      do j = 1,nch
	 vv(j,j) = vvb(j,j) - vref(j)
      enddo
      do j = 1,nch
	 do i = 1,nch
	    yy1(i,j) = ho3 * vv(i,j)
	 enddo
	 yy1(j,j) = y1(j) + yy1(j,j)
      enddo
c-----------------------------------------------------------------------
c     propagation from a to c
c-----------------------------------------------------------------------
      do j = 1,nch
	 do i = 1,nch
	    yy1(i,j) = zz(i,j) + yy1(i,j)
	 enddo
      enddo
      call invers (yy1,zz,nch,ncmx,lsym)
      do j = 1,nch
	 do i = 1,nch
	    zz(i,j) = yy4(i,j) - y3(i) * zz(i,j) * y2(j)
	 enddo
      enddo
c
c --- calculation of y1(c,b) and y4(c,b)
c
      do j = 1,nch
	 do i = 1,nch
	    yy1(i,j) = yy4(i,j)
	 enddo
      enddo
      rb = rc + hh
      call eqnsj (vvb   ,rb    ,nch   ,nblk  ,idbglm,ndmblk
     &           ,ndbmat,cr    ,enbf  ,glm   ,vqq   ,nq12
     &           ,nnq12 ,nptpot,rpot  ,vpot  ,cm    ,ndbpot
     &           ,ir1p  ,ir2p  ,idxbf ,ntrt  ,ncmx)
      call mtcopy (vv,vvb,nch,nch,ncmx,ncmx)
      do j = 1,nch
	 vv(j,j) = vvb(j,j) - vref(j)
      enddo
      do j = 1,nch
	 do i = 1,nch
	    yy4(i,j) = ho3 * vv(i,j)
	 enddo
	 yy4(j,j) = y4(j) + yy4(j,j)
      enddo
c-----------------------------------------------------------------------
c     propagation from c to b
c-----------------------------------------------------------------------
      do j = 1,nch
	 do i = 1,nch
	    yy1(i,j) = zz(i,j) + yy1(i,j)
	 enddo
      enddo
      call invers (yy1,zz,nch,ncmx,lsym)
      do j = 1,nch
	 do i = 1,nch
	    zz(i,j) = yy4(i,j) - y3(i) * zz(i,j) * y2(j)
	 enddo
      enddo
c
ccc   call matprt (zz,nch,nch,ncmx,8,1,'log-der ')
      rcur = rb
c     print *,' rcur = ',rcur
      call stepper (rcur,hcur)
      if (rcur .lt. rend-1.d-8) goto 999
      rmatch = rcur
c
      if (ibug .ge. 1) then
	 if (integrator .eq. 1) print 9001,rstart,rmatch,nstep
	 if (integrator .eq. 2) print 9002,rstart,rmatch,nstep
	 if (integrator .eq. 3) print 9003,rstart,rmatch,nstep
	 if (integrator .eq. 4) print 9004,rstart,rmatch,nstep
      endif
      if (ibug .ge. 2) then
	 call matprt (zz,nch,nch,ncmx,8,1,'log-der ')
      endif
      return
c
 9001 format ('   johnson method used between '
     &       ,f8.4,' and ',f8.4,' with ',i5,' integration steps')
 9002 format ('   johnson-manolopoulos method used between '
     &       ,f8.4,' and ',f8.4,' with ',i5,' integration steps')
 9003 format ('   johnson method used between '
     &       ,f8.4,' and ',f8.4,' with ',i5,' integration steps'
     &       ,' ; symmetric routines DSYTRx used')
 9004 format ('   johnson-manolopoulos method used between '
     &       ,f8.4,' and ',f8.4,' with ',i5,' integration steps'
     &       ,' ; symmetric routines DSYTRx used')
      end
c-----------------------------------------------------------------------
      subroutine eqnsj (hmat  ,rcur  ,nv    ,nblk  ,idbglm,ndmblk
     &                 ,ndbmat,cr    ,enbf  ,glm   ,vqq   ,nq12
     &                 ,nnq12 ,nptpot,rpot  ,vpot  ,cm    ,ndbpot
     &                 ,ir1p  ,ir2p  ,idxbf ,ntrt  ,ncmx)
      implicit double precision (a-h,o-z)
      common /steps/ step,eneff,enstb,deltae
      common /vib/ fn2q(150),fnp2q(150),fn1jj(90),fn2jj(90),fj2jj(90)
      dimension
     &     hmat(ncmx,*),idbglm(*),ndmblk(*),ndbmat(*),cr(*),enbf(*)
     &    ,glm(*),vqq(*),nptpot(*),rpot  (*),vpot  (*),cm   (*)
     &    ,ndbpot(*),ir1p(*),ir2p(*),idxbf(*)
c
      call mtinit (hmat,nv,nv,ncmx)
      call cpot (nnq12 ,nptpot,rpot  ,vpot  ,cm    ,rcur
     &          ,vqq   ,ndbpot,eneff ,enstb)
      rcur2 = rcur*rcur
      kbs = 0
      do kb = 1,nblk
	 np  = idbglm(kb)
	 ndm = ndmblk(kb)
	 ien = ndbmat(kb)
	 icr = ndbmat(kb) + nv
	 do ki = 1,ndm
	    ien = ien + 1
	    icr = icr + 1
	    do kj = 1,ki
	       sum = 0.d0
	       icon = 0
	       do iq = 1,nnq12
		  if (icon .lt. nq12) then
		     kki = idxbf(ndbmat(kb)+ki)
		     kkj = idxbf(ndbmat(kb)+kj)
		     if (((dabs(fn2jj(kki)-fn2q (iq)) .le. 0.01d0) .and.
     &                    (dabs(fn2jj(kkj)-fnp2q(iq)) .le. 0.01d0)) .or.
     &                   ((dabs(fn2jj(kkj)-fn2q (iq)) .le. 0.01d0) .and.
     &                    (dabs(fn2jj(kki)-fnp2q(iq)) .le. 0.01d0)))then
			np = np + 1
			sum = sum + glm(np) * vqq(iq)
			icon = icon + 1
		     endif
		  endif
	       enddo
	       hmat(kbs+ki,kbs+kj) = sum
	       hmat(kbs+kj,kbs+ki) = sum
	    enddo
	    hmat(kbs+ki,kbs+ki) = hmat(kbs+ki,kbs+ki)
     &                          + cr(icr) / rcur2 - enbf(ien)
	 enddo
	 kbs = kbs + ndm
      enddo
      do kt = 1,ntrt
	 ki = ir1p(kt)
	 kj = ir2p(kt)
	 hmat(ki,kj) =   cr(kt) / rcur2
	 hmat(kj,ki) =   cr(kt) / rcur2
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine stepper (rcur,hcur)
      implicit double precision (a-h,o-z)
      common /intgr/ rdb,rfn,fpt,enfdml,rmin,stbmx,fllmx,enmin,nstbmx
      common /steps/ step,eneff,enstb,deltae
c
      if (rcur .lt. rmin) then
	 step = 3.14 / fpt / dsqrt(enfdml+deltae)
      else
	 step = 3.14 / fpt / dsqrt(eneff)
      endif
      hcur = step
cc    print *,' rcur,hcur = ',rcur,hcur
      return
      end
c-----------------------------------------------------------------------
      subroutine invers (aa,bb,nn,nmx,lsym)
      implicit double precision (a-h,o-z)
      logical lsym
      parameter (ncmx1=801,lwork=64*ncmx1)
      dimension aa(nmx,*),bb(nmx,*)
      dimension ipiv(ncmx1),work(lwork)
c
      if (lsym) then
	 call dsytrf ('U',nn,aa,nmx,ipiv,work,lwork,info)
	 call dsytri ('U',nn,aa,nmx,ipiv,work,info)
	 do j = 1,nn
	    do i = j,nn
	       aa(i,j) = aa(j,i)
	    enddo
	 enddo
      else
	 call dgetrf (nn,nn,aa,nmx,ipiv,info)
	 call dgetri (nn,aa,nmx,ipiv,work,lwork,info)
      endif
      call mtcopy (bb,aa,nn,nn,nmx,nmx)
      return
      end
c-----------------------------------------------------------------------
      subroutine mtinit (aa,ii,jj,na)
      double precision aa
      dimension aa(na,*)
c
      do j = 1,jj
	 do i = 1,ii
	    aa(i,j) = 0.d0
	 enddo
      enddo
      return
c
      end
c-----------------------------------------------------------------------
      subroutine mtcopy (aa,bb,ii,jj,na,nb)
      double precision aa(na,*),bb(nb,*)
c
      do j = 1,jj
	 do i = 1,ii
	    aa(i,j) = bb(i,j)
	 enddo
      enddo
      return
c
      end
c-----------------------------------------------------------------------
      subroutine mchine (cray)
c
c     This routine is used by mtlook and matprt
c     cray = .false. is for byte machines
c     cray = .true.  is for word machines like C94 or C98
c
      logical cray
      cray = .false.
      return
c
      end
c-----------------------------------------------------------------------
      subroutine mtlook (aa,nl,nls,nc,ncs,nlmx,itype,iform,aname)
c#######################################################################
c#    quick look at a large matrix                                     #
c#    prints out a matrix of integer,real or complex numbers  a(m,n).  #
c#---------------------------------------------------------------------#
c#    aa     : matrix to be printed                                    #
c#    nl     : number of lines of the matrix                           #
c#    nls    : step for lines (lines 1,1+nls,1+2*nls, ... are printed) #
c#    nc     : number of columns of the matrix                         #
c#    ncs    : step for columns (cols 1,1+ncs,1+2*ncs,.. are printed)  #
c#    nlmx   : row dimension of aa in the calling program              #
c#    itype  : describes the type of the numbers (see matprt)          #
c#    iform  : describes the format code         (see matprt)          #
c#    name   : 8 characters to specify the name of the matrix on list. #
c#######################################################################
      logical complx,cray
      double precision aa8
      character*8      aname,anamp,blank
      dimension aa(*)
      dimension aa4(20),aa8(10),jjc(20)
      equivalence (aa4(1),aa8(1))
      data blank / '        '/
c
      if (nc .le. 0) then
	 print 9002,aname,nc
	 return
      endif
      call mchine (cray)
c
      lword = itype
      if (iform .eq. 4 .or. iform .eq. 14) go to 4000
c
      if (iform .lt. 10) ncmx = 10
      if (iform .gt. 10) ncmx =  5
      complx = itype .eq. 24 .or. itype .eq. 28
      if (complx) go to 200
      if (cray) go to 1000
      if (lword .eq. 4) go to 1000
      if (lword .eq. 8) go to 1100
c
 200  lword = itype-20
      if (cray) go to 2000
      if (lword .eq. 4) go to 2000
      if (lword .eq. 8) go to 2100
      print 9000,itype
      return
c
c --- real numbers (single precision)
c
 1000 js = -1
      ncleft = (nc-1)/ncs+1
 1010 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 1015 j = 1,ncp
	 jjc(j) = (js+j)*ncs+1
 1015 continue
      if (ncmx .eq.  5) write (6,9505) (jjc(ic),ic=1,ncp)
      if (ncmx .eq. 10) write (6,9510) (jjc(ic),ic=1,ncp)
      do 1011 i = 1,nl,nls
	 do 1012 j = 1,ncp
	    ij = i+(js+j)*nlmx*ncs
	    aa8(j) = aa(ij)
 1012    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 1011 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 1010
c
c --- real numbers (double precision)
c
 1100 js = -1
      ncleft = (nc-1)/ncs+1
 1110 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 1115 j = 1,ncp
	 jjc(j) = (js+j)*ncs+1
 1115 continue
      if (ncmx .eq.  5) write (6,9505) (jjc(ic),ic=1,ncp)
      if (ncmx .eq. 10) write (6,9510) (jjc(ic),ic=1,ncp)
      do 1111 i = 1,nl,nls
	 do 1112 j = 1,ncp
	    ij = i+(js+j)*nlmx*ncs
	    aa4(2*j-1) = aa(2*ij-1)
	    aa4(2*j)   = aa(2*ij  )
 1112    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 1111 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 1110
c
c --- complex numbers (single precision)
c
 2000 do 2010 ip = 1,2
	 if (ip .eq. 1) print 9201,aname
	 if (ip .eq. 2) print 9202,aname
	 js = -1
	 ncleft = (nc-1)/ncs+1
 2020    ncp = min0 (ncleft,ncmx)
	 if (nl .gt. 20) print 9001
	 do 2015 j = 1,ncp
	    jjc(j) = (js+j)*ncs+1
 2015    continue
	 if (ncmx .eq.  5) write (6,9505) (jjc(ic),ic=1,ncp)
	 if (ncmx .eq. 10) write (6,9510) (jjc(ic),ic=1,ncp)
	 do 2011 i = 1,nl,nls
	    do 2012 j = 1,ncp
	       ij = i+(js+j)*nlmx*ncs
	       aa8(j) = aa(2*ij+ip-2)
 2012       continue
	    call lprint (aa8,aa4,i,ncp,iform,blank)
 2011    continue
	 ncleft = ncleft-ncmx
	 js = js+ncmx
	 if (ncleft .le. 0) go to 2010
	 go to 2020
 2010 continue
      return
c
c --- complex numbers (double precision)
c
 2100 do 2110 ip = 1,2
	 if (ip .eq. 1) print 9201,aname
	 if (ip .eq. 2) print 9202,aname
	 js = -1
	 ncleft = (nc-1)/ncs+1
 2120    ncp = min0 (ncleft,ncmx)
	 if (nl .gt. 20) print 9001
	 do 2115 j = 1,ncp
	    jjc(j) = (js+j)*ncs+1
 2115    continue
	 if (ncmx .eq.  5) write (6,9505) (jjc(ic),ic=1,ncp)
	 if (ncmx .eq. 10) write (6,9510) (jjc(ic),ic=1,ncp)
	 do 2111 i = 1,nl,nls
	    do 2112 j = 1,ncp
	       ij = i+(js+j)*nlmx*ncs
	       aa4(2*j-1) = aa(2*(2*ij+ip-2)-1)
	       aa4(2*j)   = aa(2*(2*ij+ip-2)  )
 2112       continue
	    call lprint (aa8,aa4,i,ncp,iform,blank)
 2111    continue
	 ncleft = ncleft-ncmx
	 js = js+ncmx
	 if (ncleft .le. 0) go to 2110
	 go to 2120
 2110 continue
      return
c
c --- integer numbers
c
 4000 if (iform .eq.  4) ncmx = 20
      if (iform .eq. 14) ncmx = 10
      js = -1
      ncleft = (nc-1)/ncs+1
 4010 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 4015 j = 1,ncp
	 jjc(j) = (js+j)*ncs+1
 4015 continue
      if (ncmx .eq. 10) write (6,9510) (jjc(ic),ic=1,ncp)
      if (ncmx .eq. 20) write (6,9520) (jjc(ic),ic=1,ncp)
      do 4011 i = 1,nl,nls
	 do 4012 j = 1,ncp
	    ij = i+(js+j)*nlmx*ncs
	    aa4(j) = aa(ij)
 4012    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 4011 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 4010
c
 9000 format (' itype = ',i10,' should be equal to 4,8,24 or 28 ;'
     &       ,' mtlook subroutine has no action')
 9001 format ('      ')
 9002 format (' warning; mtlook did not write matrix named ',a8
     &       ,' because number of columns is ',i3)
 9010 format (a30,a8)
 9201 format (' real part of the matrix        ',a8)
 9202 format (' imaginary part of the matrix   ',a8)
 9505 format (12x,5(10x,i5,9x))
 9510 format (12x,10(4x,i5,3x))
 9520 format (12x,20(2x,i4))
      end
c-----------------------------------------------------------------------
      subroutine matprt (aa,nl,nc,nlmx,itype,iform,aname)
c#######################################################################
c#    matprt is the driver for lprint                                  #
c#    prints out a matrix of integer,real or complex numbers  a(m,n).  #
c#---------------------------------------------------------------------#
c#    aa     : matrix to be printed                                    #
c#    nl     : number of lines to be printed                           #
c#    nc     : number of columns to be printed                         #
c#    nlmx   : row dimension of aa in the calling program              #
c#    itype  : describes the type of the numbers:                      #
c#             4 for single precision or integers                      #
c#             8 for double precision or integers                      #
c#             24 (28) for single (double) precision complex           #
c#             ( a complex matrix is separated into two real matrices  #
c#               of real and imaginary parts)                          #
c#    iform  : describes the format code;                              #
c#             1 for f format code, 2 for e format code,               #
c#             3 for g format code, 4 for i format code,               #
c#             add 10 to obtain high accuracy printout                 #
c#    name   : 8 characters to specify the name of the matrix on list. #
c#######################################################################
      character*8 aname,anamp,blank
      character*30 anc(2)
      logical complx,cray
      double precision aa8
      dimension aa(*),aa4(20),aa8(10)
      equivalence (aa4(1),aa8(1))
      data blank / '        '/
     &    ,anc   / ' real part of the matrix      '
     &           , ' imaginary part of the matrix ' /
c
   10 if (nc .le. 0) then
	 print 9002,aname,nc
	 return
      endif
c
      if (nl*nc .gt. 100 000) go to 8000
      call mchine (cray)
c
      lword = itype
      if (iform .eq. 4 .or. iform .eq. 14) go to 4000
c
      if (iform .lt. 10) ncmx = 10
      if (iform .gt. 10) ncmx =  5
      complx = itype .eq. 24 .or. itype .eq. 28
      if (complx) go to 200
      if (cray) go to 1000
      if (lword .eq. 4) go to 1000
      if (lword .eq. 8) go to 1100
c
 200  lword = itype-20
      if (cray) go to 2000
      if (lword .eq. 4) go to 2000
      if (lword .eq. 8) go to 2100
      print 9000,itype
      return
c
c --- real numbers (single precision)
c
 1000 js = -1
      ncleft = nc
 1010 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 1011 i = 1,nl
	 do 1012 j = 1,ncp
	    ij = i+(js+j)*nlmx
	    aa8(j) = aa(ij)
 1012    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 1011 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 1010
c
c --- real numbers (double precision)
c
 1100 js = -1
      ncleft = nc
 1110 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 1111 i = 1,nl
	 do 1112 j = 1,ncp
	    ij = i+(js+j)*nlmx
	    aa4(2*j-1) = aa(2*ij-1)
	    aa4(2*j)   = aa(2*ij  )
 1112    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 1111 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 1110
c
c --- complex numbers (single precision)
c
 2000 do 2010 ip = 1,2
	 print 9010,anc(ip),aname
	 js = -1
	 ncleft = nc
 2020    ncp = min0 (ncleft,ncmx)
	 if (nl .gt. 20) print 9001
	 do 2011 i = 1,nl
	    do 2012 j = 1,ncp
	       ij = i+(js+j)*nlmx
	       aa8(j) = aa(2*ij+ip-2)
 2012       continue
	    call lprint (aa8,aa4,i,ncp,iform,blank)
 2011    continue
	 ncleft = ncleft-ncmx
	 js = js+ncmx
	 if (ncleft .le. 0) go to 2010
	 go to 2020
 2010 continue
      return
c
c --- complex numbers (double precision)
c
 2100 do 2110 ip = 1,2
	 print 9010,anc(ip),aname
	 js = -1
	 ncleft = nc
 2120    ncp = min0 (ncleft,ncmx)
	 if (nl .gt. 20) print 9001
	 do 2111 i = 1,nl
	    do 2112 j = 1,ncp
	       ij = i+(js+j)*nlmx
	       aa4(2*j-1) = aa(2*(2*ij+ip-2)-1)
	       aa4(2*j)   = aa(2*(2*ij+ip-2)  )
 2112       continue
	    call lprint (aa8,aa4,i,ncp,iform,blank)
 2111    continue
	 ncleft = ncleft-ncmx
	 js = js+ncmx
	 if (ncleft .le. 0) go to 2110
	 go to 2120
 2110 continue
      return
c
c --- integer numbers
c
 4000 if (iform .eq.  4) ncmx = 20
      if (iform .eq. 14) ncmx = 10
      js = -1
      ncleft = nc
 4010 ncp = min0 (ncleft,ncmx)
      anamp = aname
      if (nl .gt. 20) print 9001
      do 4011 i = 1,nl
	 do 4012 j = 1,ncp
	    ij = i+(js+j)*nlmx
	    aa4(j) = aa(ij)
 4012    continue
	 call lprint (aa8,aa4,i,ncp,iform,anamp)
	 anamp = blank
 4011 continue
      ncleft = ncleft-ncmx
      js = js+ncmx
      if (ncleft .le. 0) return
      go to 4010
 8000 print 9100,nl,nc,nlmx
      return
c
 9000 format (' itype = ',i10,' should be equal to 4,8,24 or 28 ;'
     &       ,' matprt subroutine has no action')
 9001 format ('      ')
 9002 format (' warning; matprt did not write matrix named ',a8
     &       ,' because number of columns is ',i3)
 9010 format (a30,a8)
 9100 format (' warning ---- matprt routine; nl,nc,nlmx = ',3i8
     &       ,' ; too many lines to print')
      end
c-----------------------------------------------------------------------
      subroutine lprint (aa,ia,il,nc,iform,aname)
c#######################################################################
c#    prints out a line of integer or real numbers                     #
c#---------------------------------------------------------------------#
c#    aa    : line to be printed  (real representation).               #
c#    ia    : line to be printed  (integer representation).            #
c#    nc    : number of elements in the line.                          #
c#    iform : describes the format code                                #
c#            1 for f format code, 2 for e format code,                #
c#            3 for g format code, 4 for i format code,                #
c#            add 10 to get high number of digits printout.            #
c#######################################################################
      implicit double precision (a-h,o-z)
      character*8 aname
      dimension aa(*),ia(*)
c
      igo = iform
      if (igo .gt. 10) igo = igo-6
      go to (100,200,300,400,1100,1200,1300,1400),igo
 100  print 9010, aname,il,(aa(j),j = 1,nc)
      return
 200  print 9020, aname,il,(aa(j),j = 1,nc)
      return
 300  print 9030, aname,il,(aa(j),j = 1,nc)
      return
 400  print 9040, aname,il,(ia(j),j = 1,nc)
      return
 1100 print 9110, aname,il,(aa(j),j = 1,nc)
      return
 1200 print 9120, aname,il,(aa(j),j = 1,nc)
      return
 1300 print 9130, aname,il,(aa(j),j = 1,nc)
      return
 1400 print 9140, aname,il,(ia(j),j = 1,nc)
      return
c
 9010 format (' ',a8,i4,10f12.6)
 9020 format (' ',a8,i4,1p,10d12.4)
 9030 format (' ',a8,i4,1p,10g12.5)
 9040 format (' ',a8,i4,20i6)
 9110 format (' ',a8,i4,5f24.16)
 9120 format (' ',a8,i4,1p,5d24.16)
 9130 format (' ',a8,i4,1p,5g24.16)
 9140 format (' ',a8,i4,10i12)
      end
c-----------------------------------------------------------------------
      subroutine fmmult (aa,bb,cc, ll,mm,nn, na,nb,nc)
c#######################################################################
c#    performs matrix multiplication                                   #
c#         a      = b        * c                                       #
c#          (l,n)    (l,m)      (m,n)                                  #
c#---------------------------------------------------------------------#
c#    na,nb,nc : row dimensions of a,b,c in the calling program        #
c#######################################################################
      implicit double precision (a-h,o-z)
      dimension aa(na,*),bb(nb,*),cc(nc,*)
c
      alpha = 1.d0
      beta  = 0.d0
      call dgemm ('n','n',ll,nn,mm
     &            ,alpha,bb,nb,cc,nc,beta,aa,na)
      return
c
      end
c-----------------------------------------------------------------------
      subroutine besph2 (ff,df,gg,dg,aa,arg,key,ibug)
c#######################################################################
c#    calculates linearly independent solutions of the equation        #
c#       2    2        2                                               #
c#    ( d / dx  - a / x  + 1 ) y(x) = 0                                #
c#                       -                                             #
c#    where a = l*(l+1) with l integer .ge. 0                          #
c#    + (-) sign corresponds to an open (closed) channel               #
c#    the solutions are obtained from bessel functions obtained        #
c#    in besjot, bessen, bessik subroutines                            #
c#    ff # sin (x-l*pi/2) ; gg # -cos (x-l*pi/2) for open   channels   #
c#                                                   and x infinite    #
c#---------------------------------------------------------------------#
c#    for closed channels, we have xx = dabs(arg) with arg < 0.        #
c#    then, with the definitions of abramowitz and stegun (chap. 10.2) #
c#    ff =  ( pi*x/2 )**1/2 * i     (x)                                #
c#                             l+1/2                                   #
c#    gg = -( 2*x/pi )**1/2 * k     (x)                                #
c#                             l+1/2                                   #
c#    ff #  sinh (x); gg # -exp(-x) for x  infinite                    #
c#---------------------------------------------------------------------#
c#    aa   : l*(l+1)                                                   #
c#    arg  : argument value                                            #
c#           if positive then argument z is real      (z = arg)        #
c#           if negative then argument z is imaginary (z = -i*arg)     #
c#           where i = (-1)**0.5                                       #
c#    ff,gg,df,dg : regular and irregular functions                    #
c#                  and their derivatives                              #
c#    key  : # 0 to suppress exponential factors in bessik             #
c#    ibug : .ge. 1 ,  prints function, derivative                     #
c#         : .ge. 2 ,  prints time surface integrals                   #
c#---------------------------------------------------------------------#
c#    have been tested on wronskian relation w(ff,gg) = ff*dg-df*gg= 1 #
c#    for the following values of the argument and orders :            #
c#    0 =< l < 100  and 0 =< z < 100   (error in w less than 10**-12)  #
c#    0 =< l < 30   and 0 =< z < 100*i (error in w less than 10**-10)  #
c#    outside this range check for underflows, overflows,              #
c#    divide checks and dexp capacity.                                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      dimension fff(2),dff(2)
      common /tbesp2/ time(2,2)
      data pi /3.1415 92653 58979 32384 d0/,tiny /1.d-10/
c
      if (aa .le. -0.25d0) then
	 write (6,9010) aa
	 return
      else
	 fl = -0.5d0+dsqrt (aa+0.25d0)+tiny
	 ll = fl
	 if (dabs (fl-ll) .gt. 2.d0*tiny) then
	    write (6,9010) aa
	    return
	 endif
      endif
c
      xx = dabs (arg)
      if (arg .lt. 0.d0) then
	 sk2 = -1.
	 twopim = 2.d0/pi
	 lp = ll+1
	 lm = ll-1
	 call bessik (ll,xx,bi ,bk ,key)
	 call bessik (lm,xx,bim,bkm,key)
	 call bessik (lp,xx,bip,bkp,key)
	 dbi =  (ll*bim+lp*bip)/(ll+lp)
	 dbk = -(ll*bkm+lp*bkp)/(ll+lp)
	 ff =   xx*bi
	 df =  (xx*dbi+bi)
	 gg = -twopim*xx*bk
	 dg = -twopim*(xx*dbk+bk)
      else
	 sk2 = 1.
	 call besjot (ll,xx,bj,dbj,r)
	 call bessen (ll,xx,bn,dbn,rg)
	 do 210 i = 1,ll
	    bj = bj/r
	    dbj = dbj/r
	    bn = bn*rg
	    dbn = dbn*rg
 210     continue
	 ff = xx*bj
	 df = xx*dbj+bj
	 gg = xx*bn
	 dg = xx*dbn+bn
      endif
c
c --- time surface integrals calculation
c
      fff(1) = ff
      fff(2) = gg
      dff(1) = df
      dff(2) = dg
      vv = aa/(xx*xx)-sk2
      do 1010 i = 1,2
	 do 1011 j = 1,2
	    time(i,j) = -0.125*(fff(i)*dff(j)+dff(i)*fff(j)
     &                  +2.*xx*(fff(i)*vv*fff(j)-dff(i)*dff(j)) )
 1011    continue
 1010 continue
c
      if (ibug .gt. 0) then
	 wm1 = ff*dg-df*gg-1.d0
	 write (6,9000) aa,arg,ff,df,gg,dg,wm1
	 if (ibug .gt. 1) then
	    call matprt (time,1,4,1,8,12,'time    ')
	 endif
      endif
      return
c
 9000 format (' besph2 ',1f12.4,1f12.4,1p,4d20.12,1d9.1)
 9010 format (' ****** besph2 routine; aa = ',f12.2,' is not l*(l+1)'
     &       ,' with l integer .ge. 0; request aborted ******')
      end
c-----------------------------------------------------------------------
      subroutine bessik (ll,xx,bi,bk,key)
c#######################################################################
c#    calculation of modified spherical bessel functions of the third  #
c#    kind  :                                                          #
c#    (pi/(2*x))**1/2 * i     (x)  by formula 10.2.5 (p. 443)          #
c#                       l+1/2                                         #
c#    (pi/(2*x))**1/2 * k     (x)  by formula 10.2.15 (p. 444)         #
c#                       l+1/2                                         #
c#    (definitions of abramowitz and stegun, chap. 10)                 #
c#---------------------------------------------------------------------#
c#    ll    : order (integer)                                          #
c#    xx    : argument (positive real number)                          #
c#    bi,bk : bessel functions  ; key # 0 to suppress exp. factors     #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /cfaclog/ fct(nfctmx)
      data pi /3.1415 92653 58979 32384 d0/ ,ntimes /0/
c
      ntimes  = ntimes+1
      if (ntimes .eq. 1) call faclog
      bi = 0.d0
      bk = 0.d0
      if (ll .lt. 0) return
c
      expx = dexp (-xx)
c
      fl = ll
      flog2 = dlog (2.d0)
      flogx = dlog (xx)
      ci = dexp (fl*flogx)
      floga = dlog (0.5d0)+2.d0*flogx
      sum = 0.d0
      kk = 0
      lk = ll
 100  st = kk*floga-fct(kk+1)-( fct(2*lk+2)-lk*flog2-fct(lk+1) )
      kk = kk+1
      lk = ll+kk
      term = dexp(st)
      sum = sum+term
      if ((term/sum) .gt. 1.d-12) go to 100
      bi = ci*sum
c
      sum = 0.5d0/xx
      if (ll .eq. 0) go to 200
	 shox = dlog (0.5d0/xx)
	 shox1 = 0.d0
	 sum = 0.d0
	 do 201 kk = 0,ll
	    shox1 = shox+shox1
	    st = fct(ll+kk+1)-fct(kk+1)-fct(ll-kk+1)+shox1
	    if (st .lt. -180.d0) go to 200
c           if (st .lt. -500.d0) go to 200
	    sum = sum+dexp(st)
 201     continue
 200  bk = pi*expx*sum
c
      if (key .eq. 0) return
	 bi = bi*expx
	 bk = bk/expx
      return
c
      end
c-----------------------------------------------------------------------
      subroutine besjot (l,x,f,df,r)
c
c   this program generates the standard and modified versions of the sph
c   bessel-functions of first(besjot)- and second(bessen)-kind respectiv
c   for definitions compare: 'nbs-h1ndbook of mathematical functions',
c   (abramowitz+stegun,eds./n.y.:1964), sections 10.1.1 on page 437 for
c   standard versions and ss.10.2.2 + 10.2.3 on p.443 for the modified o
c   l=index(natural numbers including zero), x=argument(real,d.p.), f=ou
c
c   the sign of the argument is used to determine the versions:
c   the outcomes f(=first-kind-functions) and g(=second-kind-f.) must be
c   divided (resp. multiplied) by the l-th power of the reduction logfac
c   to get the modified versions, use argument with negative sign }
c
c   by formulas 10.1.31 on page 439 loc.cit. and 10.2.7 on p.443 ibid.,
c   solutions have been tested to be correct to twelve places at least i
c   range combining x=1...441 and l=0...340 .
c
c   besjot is divided into three parts, corresponding to wether x > l, o
c   while x < l, beeing 0.5*x*x < 2*l or 0.5*x*x > 2*l respectively .
c
c
c modif. pour x plus grand que l dans besjot  : r=1
c        qq soit x dans bessen : r=1
c   pour eviter les overflows ou underflows dans le prog. appele pour 50
c      version jan. 77
c
      implicit double precision (a-h,o-z)
c
       f = 0.d0
       r = 1.d0
      if(x)  51,50,52
 50   if(l.eq.0)  f=1.d0
      df = 0.d0
      return
c
 51   sinix = dsinh(-x)
      cosix = dcosh(-x)
       w = -1.d0
      goto 53
 52   pi = 6.283185307179586d0
      xr = dmod(x,pi)
      sinix = dsin(xr)
      cosix = dcos(xr)
       w = +1.d0
 53    z = 1.d0 / dabs(x)
       a = dfloat(l)
       r = a * z
      if(dabs(x)-a) 2,2,1
   2  if(0.5d0*x*x-2.d0*a) 4,4,3
c
c   for the following version see page 439, section 10.1.19 loc.cit., an
c   section 10.1.11 on page 438 ibidem}
c   for the modified case look up ss.10.2.18 and 10.2.13 on ps.444 and 4
c   this version is used, if x > l .
c
  1   f0 = sinix * z
      g0 =-cosix * z
	r=1.d0
      if(l) 11,11,12
 11    f = f0
      df =-g0 - f0*z
      return
 12   if(l-1) 13,13,14
 13    f = w * (f0-cosix) * z
      df = f0   - 2.0d0*f*z
      return
 14   f1 = w * (f0-cosix) * z
      if(l.eq.2)  goto 15
       j = l-2
      do 10 i=1,j
      f2 = w * (f1*dfloat(2*i+1)*z - f0  )
      f0 = f1
   10 f1 = f2
   15  f = w * (f1*dfloat(2*l-1)*z - f0  )
      df = f1   - (l+1)*f*z
      return
c
c   for the following version see page 453, example 2 loc.cit.}
c   this version is used, if x < l and 0.5*x*x > 2*l .
c
  3    n = a + 25.d0 + dsqrt(a)
      b0 = 0.d0
      b1 = 1.d0
      do 20 j=1,n
      b2 = w * (b1*dfloat(2*(n-j)+3)*z - b0/r) / r
      b0 = b1
      if(n-l-j) 22,21,22
 21    f = b2
      goto 20
 22   if(n-l+1-j)  20,23,20
 23   df = b2
   20 b1 = b2
      df = w**(l-1) * (df/b1) * sinix*z
       f = w**l * (f/b1) * sinix*z
      df = df*r - (l+1)*f*z
      return
c
c   for the following version see p1ge 437, formula 10.1.2  loc.cit.}
c   for the modified case formula 10.2.5 on page 443 is valid.
c   this version is used, if x < l and 0.5*x*x < 2*l .
c
c
  4    y = -w * 0.5d0 * x * x
      s0 = dfloat(2*l-1)
      s1 = dfloat(2*l+1)
      p0 = 1.d0
      p1 = 1.d0
      c0 = 1.d0
      c1 = 1.d0
      do 30 i=1,15
      s0 = s0 + 2.d0
      s1 = s1 + 2.d0
      p0 = y*p0/(s0*dfloat(i))
      p1 = y*p1/(s1*dfloat(i))
      c0 = c0 + p0
   30 c1 = c1 + p1
       q = 1.d0
      if(l.eq.1)  goto 32
       j = l - 1
      do 31 i=1,j
   31  q = q * a/dfloat(2*i+1)
       f = q * a/dfloat(2*l+1) * c1
      df = q*c0*r - dfloat(l+1)*f*z
      return
  32   f = c1 / 3.d0
      df = c0*r - 2.d0*f*z
      return
c
c
	   entry bessen(l,x,g,dg,r)
c
c   spherical bessel-(and modified bessel-) functions of the second kind
c   this version is valid for all indices and arguments.
c
       g = 0.d0
       a = dfloat(l)
       r = 1.d0
      if(x) 61,60,62
 60   write(6,300)
300   format(1h0,'*******  argument of spherical bessel-function of seco
     1nd kind should not be zero }')
      return
 61   sinix = dsinh(-x)
      cosix = dcosh(-x)
       w = -1.d0
      goto 63
 62   pi = 6.283185307179586d0
      xr = dmod(x,pi)
      sinix = dsin(xr)
      cosix = dcos(xr)
       w = +1.d0
 63    z = 1.d0 / dabs(x)
      g0 =-w * cosix * z
      f0 = w * sinix * z
      if(l) 41,41,42
 41    g = g0
      dg = w*f0 - g0*z
      return
 42        r=1.d0
       if(l-1) 43,43,44
 43    g = w * (g0-sinix)*z
      dg = g0   - 2.0d0*g*z
      return
 44   g1 = w * (g0-sinix)*z
      if(l.eq.2)  goto 45
       j = l-2
      do 40 i=1,j
      g2 = w * (g1*dfloat(2*i+1)*z - g0  )
      g0 = g1
   40 g1 = g2
  45   g = w * (g1*dfloat(2*l-1)*z - g0  )
      dg = g1   - dfloat(l+1)*g*z
      return
      end
c-----------------------------------------------------------------------
      double precision function f3j (fj1,fj2,fj3, fm1,fm2,fm3)
c#######################################################################
c#    calculates 3j coefficients from racah formula                    #
c#    (messiah: t2, p 910; formula 21) .                               #
c#    clebsch-gordan coefficients are given by (p. 908, formula 12) :  #
c#                         j -j +m                |j    j     j|       #
c#    <j j m m |j m> = (-1) 1  2   (2*j+1)**(0.5) | 1    2     |       #
c#      1 2 1 2                                   |m    m    -m|       #
c#                                                | 1    2     |       #
c#######################################################################
      implicit double precision (a-h,o-z)
      integer t,tmin,tmax
      parameter (nfctmx=1001)
      data tiny,zero,one /0.01d0,0.d0,1.d0/ ,ntimes /1/
      common /cfaclog/ fct(nfctmx)
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      cc = zero
      if (fj3 . gt. (fj1+fj2+tiny))      go to 100
      if (dabs(fj1-fj2) .gt. (fj3+tiny)) go to 100
      if (dabs(fm1+fm2+fm3) .gt. tiny)   go to 100
      if (dabs(fm1) .gt. (fj1+tiny))     go to 100
      if (dabs(fm2) .gt. (fj2+tiny))     go to 100
      if (dabs(fm3) .gt. (fj3+tiny))     go to 100
      fk1 = fj3-fj2+fm1
      fk2 = fj3-fj1-fm2
      fk3 = fj1-fm1
      fk4 = fj2+fm2
      fk5 = fj1+fj2-fj3
      fk1m = fk1-tiny
      fk2m = fk2-tiny
      fk1p = fk1+tiny
      fk2p = fk2+tiny
      if (fk1m .lt. zero) k1 = fk1m
      if (fk1p .gt. zero) k1 = fk1p
      if (fk2m .lt. zero) k2 = fk2m
      if (fk2p .gt. zero) k2 = fk2p
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = fj1+fj2-fj3+one+tiny
      n2 = fj2+fj3-fj1+one+tiny
      n3 = fj3+fj1-fj2+one+tiny
      n4 = fj1+fm1+one+tiny
      n5 = fj2+fm2+one+tiny
      n6 = fj3+fm3+one+tiny
      n7 = fj1-fm1+one+tiny
      n8 = fj2-fm2+one+tiny
      n9 = fj3-fm3+one+tiny
      n10 = fj1+fj2+fj3+2.d0+tiny
      x = fct(n1)+fct(n2)+fct(n3)+fct(n4)+fct(n5)+fct(n6)
     &   +fct(n7)+fct(n8)+fct(n9)-fct(n10)
      x = 0.5d0*x
      do 10  t = tmin,tmax
	 phase = one
	 if (mod(t,2) .ne. 0) phase = -one
	 cc = cc+phase*dexp(-fct(t+1)   -fct(k1+t+1)-fct(k2+t+1)
     &                      -fct(k3-t+1)-fct(k4-t+1)-fct(k5-t+1)+x)
 10   continue
      fsp = dabs(fj1-fj2-fm3)+tiny
      ns = fsp
      if (mod(ns,2) .gt. 0) cc = -cc
 100  f3j = cc
      return
      end
c-----------------------------------------------------------------------
      double precision function f6j (fj1,fj2,fj3,fl1,fl2,fl3)
c#######################################################################
c#    calculation of 6j-coefficients                                   #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /cfaclog/ fct(nfctmx)
      data tiny /.01/ ,ntimes /1/
c
      if (ntimes .eq. 1) call faclog
      ntimes = ntimes+1
      d = fdelta (fj1,fj2,fj3)
      d = d*fdelta (fj1,fl2,fl3)
      d = d*fdelta (fl1,fj2,fl3)
      d = d*fdelta (fl1,fl2,fj3)
      f6j = 0.d0
      if (dabs(d) .eq. 0.d0) return
c
      fk1 = fj1+fj2+fj3
      fk2 = fj1+fl2+fl3
      fk3 = fl1+fj2+fl3
      fk4 = fl1+fl2+fj3
      fk5 = fj1+fj2+fl1+fl2
      fk6 = fj2+fj3+fl2+fl3
      fk7 = fj3+fj1+fl3+fl1
      fmin = dmin1 (fk5,fk6,fk7)
      fmax = dmax1 (fk1,fk2,fk3,fk4)
      min = fmin+tiny
      max = fmax+tiny
      k1 = fk1+tiny
      k2 = fk2+tiny
      k3 = fk3+tiny
      k4 = fk4+tiny
      k5 = fk5+tiny
      k6 = fk6+tiny
      k7 = fk7+tiny
      if (min-max) 1000,3,3
 3    if (max) 1000,4,4
 4    if (min) 1000,5,90
 5    k1 = -k1
      k2 = -k2
      k3 = -k3
      k4 = -k4
      bot = fct(k1+1)+fct(k2+1)+fct(k3+1)+fct(k4+1)+fct(k5+1)+fct(k6+1)
     &     +fct(k7+1)
      bot = dexp(bot)
      f6j = d/bot
      return
c
 90   f6j = 0.
      do 100 i = max,min
	 boite = ((-1.)**i)
	 iz = i+1
	 m1 = i-k1
	 m2 = i-k2
	 m3 = i-k3
	 m4 = i-k4
	 m5 = k5-i
	 m6 = k6-i
	 m7 = k7-i
	 dot = fct(iz+1)
	 bot = fct(m1+1)+fct(m2+1)+fct(m3+1)+fct(m4+1)+fct(m5+1)
     &        +fct(m6+1)+fct(m7+1)
	 b1 = dot-bot
	 boite = boite*dexp(b1)
	 f6j = f6j+boite
 100  continue
      f6j = f6j*d
 1000 return
c
      end
c-----------------------------------------------------------------------
      double precision function fdelta (fl1,fl2,fl3)
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /cfaclog/ fct(nfctmx)
      data   eps /.01/
c
      ia=fl1+fl2+fl3 +eps
      a=2.*(fl1+fl2+fl3)+1.
      ib=a +eps
      ib=ib/2
      if(ib-ia)1,6,1
    6 continue
      ik1=fl1+fl2-fl3+eps
      ik2=fl2+fl3-fl1+eps
      ik3=fl3+fl1-fl2+eps
      kk=fl1+fl2+fl3+1+eps
      if(ik1)1,2,2
    2 if(ik2)1,3,3
    3 if(ik3)1,4,4
    4 d1=fct(kk+1)
      d2=fct(ik1+1)+fct(ik2+1)+fct(ik3+1)
      d3 = (d2 - d1) / 2.d0
      fdelta = dexp (d3)
      go to 5
    1 fdelta=0.
    5 return
c
      end
c-----------------------------------------------------------------------
      double precision function f9j (fj1   ,fj2   ,fj12
     &                              ,fj3   ,fj4   ,fj34
     &                              ,fj13  ,fj24  ,fj)
c#######################################################################
c#    calculates 9j coefficients from 6j coefficients.                 #
c#    formula (41) p.917 messiah (t.2)                                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      data tiny,tiny1 /1.d-2,1.d-40/
      g1 = dabs (fj1-fj)
      g2 = dabs (fj2-fj34)
      g3 = dabs (fj3-fj24)
      g4 = fj1+fj
      g5 = fj2+fj34
      g6 = fj3+fj24
      gmn = dmax1 (g1,g2,g3)
      gmx = dmin1 (g4,g5,g6)
      sum = 0.d0
      if (gmn .gt. (gmx+tiny)) go to 10
      g = gmn
 20   a1 = f6j (fj1   ,fj2   ,fj12  ,fj34  ,fj   ,g)
      if (dabs (a1) .lt. tiny1) go to 30
      a2 = f6j (fj3   ,fj4   ,fj34  ,fj2   ,g    ,fj24)
      if (dabs (a2) .lt. tiny1) go to 30
      a3 = f6j (fj13  ,fj24  ,fj    ,g     ,fj1  ,fj3)
      if (dabs (a3) .lt. tiny1) go to 30
      n = 2*g+tiny
      phase = 1.d0
      if (mod(n,2) .gt. 0) phase = -1.d0
      sum = sum+phase*(2*g+1)*a1*a2*a3
 30   g = g+1
      if (g .lt. (gmx+tiny)) go to 20
 10   f9j = sum
      return
      end
c-----------------------------------------------------------------------
      subroutine faclog
c#######################################################################
c#    initialisation of logarithms of factorials array                 #
c#######################################################################
      implicit double precision (a-h,o-z)
      parameter (nfctmx=1001)
      common /cfaclog/ fct(nfctmx)
      data ntimes /0/
c
      ntimes = ntimes+1
      if (ntimes .gt. 1) return
      fct(1) = 0.d0
      do 10 i = 1,nfctmx-1
	 ai = i
	 fct(i+1) = fct(i)+dlog(ai)
 10   continue
c
      return
      end
c-----------------------------------------------------------------------
      include 'cpot.f'
c-----------------------------------------------------------------------

