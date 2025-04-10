!!!*************************************************************
! 文件/File: Ay.gensub.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: Ay.gensub.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

      module gensub
      contains
********************************************************************************
C      gensub.f  -- standard subroutines, to calculate wigner coeffs,
c			interpolation, root-finding, numerical integration
c			Coulomb function programs also.
C**********************************************************C
C
C WIGNER 3J, 6J, AND 9J FUNCTION SUBPROGRAMS.
C NOTE:  ALL ANGULAR MOMENTUM QUANTUM NUMBERS SUPPLIED TO THESE
C        FUNCTIONS ARE INTEGERS WHICH ARE TWICE THE VALUE OF THE
C        ACTUAL ANGULAR MOMENTUM DESIRED.  (THIS ALLOWS FOR HALF-
C        INTEGER VALUES CONVENIENTLY.)  ALSO YOU MUST CALL SETUP
C        ONE TIME BEFORE CALLING ANY OF THESE SO THAT THE RELE-
C        VANT FACTORIALS CAN BE CALCULATED ONCE AND FOR ALL AND
C        STORED, AS IN THE ABOVE EXAMPLE.  
C
C**********************************************************C
C                                                                       
      function coeff(l1,l2,l1p,l2p,l,k)
      implicit real*8(a-h,o-z)
!      implicit none 
!      real*8 :: front,coeff,t1,t2,s1,thrj,sixj
!      integer :: l1,l2,l1p,l2p,iz,k,kd,l1d,l1pd,l2d,l2pd,l,ld
!      external :: thrj,sixj
      front=(2*l1+1)*(2*l2+1)*(2*l1p+1)*(2*l2p+1)
      front=dsqrt(front)*(-1)**(l1+l1p+l)
      l1d=2*l1
      l2d=2*l2
      l1pd=2*l1p
      l2pd=2*l2p
      ld=2*l
      kd=2*k
      iz=0
      t1=thrj(l1d,kd,l1pd,iz,iz,iz)
      t2=thrj(l2d,kd,l2pd,iz,iz,iz)
      s1=sixj(l1d,l2d,ld,l2pd,l1pd,kd)
      coeff=front*t1*t2*s1
      end function coeff
c
      FUNCTION XNINEJ(J11,J12,J13,J21,J22,J23,J31,J32,J33)              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      KMIN1 = IABS(J11-J33)                                             
      KMIN2 = IABS(J32-J21)                                             
      KMIN3 = IABS(J23-J12)                                             
      KMAX1 = J11+J33                                                   
      KMAX2 = J32+J21                                                   
      KMAX3 = J23+J12                                                   
C                                                                       
      IF(KMIN2.GT.KMIN1) KMIN1=KMIN2                                    
       IF(KMIN3.GT.KMIN1) KMIN1=KMIN3                                   
      IF(KMAX2.LT.KMAX1) KMAX1=KMAX2                                    
      IF(KMAX3.LT.KMAX1) KMAX1=KMAX3                                    
C                                                                       
      KMIN1 = KMIN1 + 1                                                 
      KMAX1 = KMAX1 + 1                                                 
       XNINEJ = 0.D0                                                    
      IF(KMIN1.GT.KMAX1) GO TO 1000                                     
      DO 100 K1 = KMIN1,KMAX1,2                                         
      K = K1 - 1                                                        
      S1 = SIXJ(J11,J21,J31,J32,J33,K)                                  
      S2 = SIXJ(J12,J22,J32,J21,K,J23)                                  
      S3 = SIXJ(J13,J23,J33,K,J11,J12)                                  
      P = (K+1)*((-1)**K)                                               
      XNINEJ = XNINEJ + P*S1*S2*S3                                      
  100 CONTINUE                                                          
 1000 CONTINUE                                                          
      RETURN                                                            
      END  function xninej
C
      FUNCTION THRJ(J1D,J2D,J3D,M1D,M2D,M3D)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      X1 = J1D/2.D0                                                     
      X2 = J2D/2.D0                                                     
      X3 = J3D/2.D0                                                     
      Y1 = M1D/2.D0                                                     
      Y2 = M2D/2.D0                                                     
      Y3 = M3D/2.D0                                                     
C                                                                       
C -- NEXT COME THE TRIANGULARITY CHECKS:
C
      IF(J1D+J2D-J3D.LT.0) GO TO 9998                                   
      IF(J2D+J3D-J1D.LT.0) GO TO 9998                                   
      IF(J3D+J1D-J2D.LT.0) GO TO 9998                                   
      IF(J3D.LT.IABS(J1D-J2D)) GO TO 9998
      IF(M1D+M2D+M3D.NE.0) GO TO 9998
      LLL = J1D+J2D+J3D
      IF(2*(LLL/2) .NE. LLL) GO TO 9998
C
      KMIN = (J3D-J1D-M2D)/2                                            
      KMIN1 = KMIN                                                      
      KMIN2 = (J3D-J2D+M1D)/2                                           
      IF(KMIN2.LT.KMIN) KMIN=KMIN2                                      
      KMIN = (-1)*KMIN                                                  
      KMAX = X1+X2-X3 +0.1D0                                            
      KMAX1 = KMAX                                                      
      KMAX2 = X1-Y1                                                     
      KMAX3 = X2+Y2                                                     
      IF(KMAX2.LT.KMAX) KMAX=KMAX2                                      
      IF(KMAX3.LT.KMAX) KMAX=KMAX3                                      
      IF(KMIN.LT.0) KMIN = 0                                            
      IF(KMIN.GT.KMAX) GO TO 9998                                       
C                                                                       
      JMIN = KMIN+1                                                     
      JMAX = KMAX+1                                                     
      TERM1 = FRONTL(X1,X2,X3,Y1,Y2,Y3)                                 
	iphase=iabs((j1d-j2d-m3d)/2)
	msign=(-1)**iphase
cg     MSIGN = (-1)**((J1D-J2D-M3D)/2)                                   
      SUM = 0.D0                                                        
      DO 10 I1 = JMIN,JMAX                                              
      I = I1 - 1                                                        
      TERM2 = FL(I1)+FL(KMIN1+I1)+FL(KMIN2+I1)                          
      TERM2 = TERM2+FL(KMAX1-I+1)+FL(KMAX2-I+1)+FL(KMAX3-I+1)           
      TERM= DEXP(TERM1-TERM2)                                           
      TERM = TERM*MSIGN*((-1)**I)                                       
      SUM = SUM + TERM                                                  
  10  CONTINUE                                                          
      THRJ = SUM                                                        
      GO TO 9999                                                        
 9998 THRJ = 0.D0                                                       
 9999 CONTINUE                                                          
      RETURN                                                            
      END   function thrj
C                                                                       
      FUNCTION FL(I)                                                    
       IMPLICIT REAL*8(A-H,O-Z)                                         
CCC   DIMENSION FACL(60)                                                
      COMMON/FACTOR/FACL(200)                                           
      FL = FACL(I)                                                      
      RETURN                                                            
      END function fl 
C
C** ** **                                                               
C-- THIS SUBROUTINE INITIALIZES BY FINDING THE LOGARITHM
C---OF THE FIRST 199 FACTORIALS AND STORING THEM.
C
      SUBROUTINE SETUP                                                  
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/FACTOR/FACL(200)                                           
      N = 199                                                           
      FACL(1) = 0.D0                                                    
      DO 100 I = 2,N                                                    
      I1 = I - 1                                                        
      FACL(I) = FACL(I1) + DLOG(I1*1.D0)                                
  100 CONTINUE                                                          
      RETURN                                                            
      END  subroutine setup
C** ** **                                                               
C                                                                       
      FUNCTION FRONTL(X1,X2,X3,Y1,Y2,Y3)                                
      IMPLICIT REAL*8(A-H,O-Z)                                          
      L1 = X1+X2-X3 +1.1D0                                              
      L2 = X2+X3-X1 +1.1D0                                              
      L3 = X3+X1-X2 +1.1D0                                              
      L4 = X1+X2+X3+1+1.1D0                                             
      L5 = X1+Y1+1.1D0                                                  
      L6 = X1-Y1+1.1D0                                                  
      L7 = X2+Y2+1.1D0                                                  
      L8 = X2-Y2+1.1D0                                                  
      L9 = X3+Y3+1.1D0                                                  
      L10 = X3-Y3+1.1D0                                                 
      FRONTL = FL(L1)+FL(L2)+FL(L3)-FL(L4)+FL(L5)+FL(L6)                
      FRONTL = FRONTL +FL(L7)+FL(L8)+FL(L9)+FL(L10)                     
      FRONTL = FRONTL/2.D0                                              
      RETURN                                                            
      END  function frontl
C
      FUNCTION SIXJ(J1D,J2D,J3D,J4D,J5D,J6D)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
C                                                                       
C -- CHECK THAT TRIANGULARITY CONDITIONS ARE SATISFIED.
C
      IF(J3D.LT.IABS(J1D-J2D) .OR. J3D.GT.J1D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J4D-J2D) .OR. J6D.GT.J4D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J1D-J5D) .OR. J6D.GT.J1D+J5D) GO TO 9998
      IF(J3D.LT.IABS(J4D-J5D) .OR. J3D.GT.J4D+J5D) GO TO 9998
      K1=J1D+J2D+J3D
      K2=J4D+J2D+J6D
      K3=J6D+J1D+J5D
      K4=J3D+J4D+J5D
      IF(2*(K1/2).NE.K1 .OR. 2*(K2/2).NE.K2) GO TO 9998
      IF(2*(K3/2).NE.K3 .OR. 2*(K4/2).NE.K4) GO TO 9998
C
C -- NOW GO AHEAD AND CALCULATE THE SIXJ.
C
      JM1 = (J1D+J2D+J3D)/2                                             
      JM2 = (J1D+J5D+J6D)/2                                             
      JM3 = (J4D+J2D+J6D)/2                                             
      JM4 = (J4D+J5D+J3D)/2                                             
      JX1 = (J1D+J2D+J4D+J5D)/2                                         
      JX2 = (J2D+J3D+J5D+J6D)/2                                         
      JX3 = (J3D+J1D+J6D+J4D)/2                                         
C                                                                       
      JM = JM1                                                          
      IF(JM2.GT.JM) JM = JM2                                            
      IF(JM3.GT.JM) JM = JM3                                            
      IF(JM4.GT.JM) JM = JM4                                            
      JX = JX1                                                          
      IF(JX2.LT.JX) JX = JX2                                            
      IF(JX3.LT.JX) JX = JX3                                            
      KM = JM+1                                                         
      KX = JX+1                                                         
      IF(KM.GT.KX) GO TO 9998                                           
      TERM1 = FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)                           
      SIXJ = 0.D0                                                       
      DO 10 I1 = KM,KX                                                  
      I = I1 - 1                                                        
      TERM2 = FL(I+2)-FL(I+1-JM1)-FL(I+1-JM2)-FL(I+1-JM3)               
      TERM2 = TERM2-FL(I+1-JM4)-FL(JX1-I+1)-FL(JX2-I+1)                 
      TERM2 = TERM2 -FL(JX3-I+1)                                        
      TERM = DEXP(TERM1+TERM2) * ((-1)**I)                              
      SIXJ = SIXJ + TERM                                                
   10 CONTINUE                                                          
      GO TO 9999                                                        
 9998 CONTINUE                                                          
      SIXJ = 0.D0                                                       
 9999 CONTINUE                                                          
      RETURN                                                            
      END  function sixj
C                                                                       
      FUNCTION FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      FRTSJL = DL(J1D,J2D,J3D) + DL(J1D,J5D,J6D)                        
      FRTSJL = FRTSJL + DL(J4D,J2D,J6D) + DL(J4D,J5D,J3D)               
      RETURN                                                            
      END function frtsjl
C                                                                       
      FUNCTION DL(J1D,J2D,J3D)                                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      L1 = (J1D+J2D-J3D)/2                                              
      L2 = (J2D+J3D-J1D)/2                                              
      L3 = (J3D+J1D-J2D)/2                                              
      L4 = (J1D+J2D+J3D)/2 + 1                                          
      DL = FL(L1+1)+FL(L2+1)+FL(L3+1)-FL(L4+1)                          
      DL = DL/2.D0                                                      
      RETURN                                                            
      END  function dl 
C
      FUNCTION CLEBSH(J1D,J2D,JD,M1D,M2D,MD)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      CLEBSH = 0.D0                                                     
      IF(M1D+M2D.NE.MD) GO TO 100                                       
      MMD = -MD                                                         
      Q = THRJ(J1D,J2D,JD,M1D,M2D,MMD)                                  
! Modified by ST 16/12/04
      PHASE = ((-1)**(20 + (J2D-J1D-MD)/2))                             
      CLEBSH = Q*PHASE*DSQRT(JD+1.D0)                                   
!      phase= (-1)**((j1d-j2d+md)/2)
!      clebsh=q*phase*sqrt(jd+1.d0)
  100 CONTINUE                                                          
      RETURN                                                           
      END  function clebsh
************************************************************************
***********************************************************************
c*****************************************************************************
        double precision function delt(i,j)
        implicit double precision (a-h,o-z)
        delt=0.d0
        if(i.eq.j) delt=1.d0
        return
        end function delt
c***********************************************************************
        subroutine root(xtry,hh,iroute,ep,xest,yest,yfcn)
        implicit real*8(a-h,o-z)


c***  subroutine root locates a root of yfcn(x,y=0),
c***  returning finally (xest,yest)
c***  xtry is initial guess,hh is step size for the iteration.
c***  iroute=0 means hh is kept constant in the initial (coarse)
c***  root location loop.
c***  iroute=1 means hh can change sign but not magnitude.
c***  iroute=2 means both magnitude & sign of hh can vary.
c***  after the initial loop, a 3-point method is used to
c***  speed convergence.
c***
        iq = 3
        u1 = xtry
	f1=yfcn(u1)
        mtest = 0
        do 420 ki = 1,30
        if(ki.eq.2 .and. iroute.eq.2) go to 421
        if(mtest.ne.0) go to 920
        u2 = u1 + hh
  421   f2=yfcn(u2)
        if(f1*f2 .le. 0.d0) mtest = 1000
        if(mtest.ne.0)  go to 3055
        q1 = dabs(f1)
        q2 = dabs(f2)
        if(ki.eq.1 .and. iroute.eq.2) go to 3047
        if(iroute-1) 3033,3044,3044
 3044   if(q1.ge.q2 .and. mtest.eq.0) u1=u2
        if(q1.ge.q2 .and. mtest.eq.0) f1=f2
        if(q2.gt.q1 .and. mtest.eq.0) hh = (-1)*hh
        go to 3055
 3033   if(mtest.eq.0) u1 = u2
        if(mtest.eq.0) f1 = f2
        go to 3055
 3047   slp = (f2-f1)/(u2-u1)
        uest = u1 - f1/slp
        if(q2.lt.q1) f1 = f2
        if(q2.lt.q1) u1 = u2
        u2 = uest
 3055   continue
  420   continue
        if(mtest.ne.1000)  go to 8500
c
  920   continue
ccccccccccccccccccccccccccccccccc      ep = 1.d-06
        mchk = 0
        do 440 i = 1,50
        if(mchk.ne.0) go to 8600
        if(i.gt.iq) go to 429
        slp = (f2-f1)/(u2-u1)
        uest = u1 - f1/slp
        if(i.lt.iq) uest = (u1+u2)/2.d0
        u3 = uest
        q3 = u3
        q2 = u2
        q1 = u1
  429   fest=yfcn(uest)
        f3 = fest
        u3 = uest
        if(i.ge.iq) go to 439
        if(fest*f1 .lt. 0.d0) f2=fest
        if(fest*f1 .lt. 0.d0) u2=uest
        if(fest*f1 .ge. 0.d0) f1=fest
        if(fest*f1 .ge. 0.d0) u1=uest
        go to 440
  439   continue
c
        aa = (u2-u1)*f3+(u1-u3)*f2+(u3-u2)*f1
        denom = (u3-u2)*(u2-u1)*(u2-u1)
        aa = aa/denom
        bb = (u2-u1)*(2.d0*u3-u2-u1)*f3
        bb = bb- (u3-u1)*(u3-u1)*f2+(u3-u2)*(u3-u2)*f1
        bb = bb/denom
        cc = (u3-u1)*f3/(u2-u1)
        q0q = dabs(bb*bb - 4.d0*aa*cc)
        q0q = dsqrt(q0q)
        den1 = bb+q0q
        den2 = bb-q0q
        qden1 = dabs(den1)
        qden2 = dabs(den2)
        if(qden2.ge.qden1) den1 = den2
        u1 = u2
        u2 = u3
        u3 = u2 - 2.d0*cc/den1
        f1 = f2
        f2 = f3
        tst = dabs(u3-u2)
        if(tst.lt.ep) mchk = 1000
        uest = u3
  440   continue
        if(mchk.ne.1000) go to 8500
        go to 8600
 8500   write(6,8501)
 8501   format(3x,'no convergence was achieved')
 8600   continue
        xest = uest
        yest = fest
        return
        end subroutine root
***********************************************************************
c ********
c	C. Greene, 6-8-87   -- modified 11-9-87
c  The following is a main program which shows how to use subroutine seaton.
c
c       The form of the call is:
c
c	     call seaton(l,eryd,r,zion,f,fp,g,gp)
c
c	Here l=orbital angular momentum
c	     eryd = electron energy in rydbergs
c            r=radial distance
c 	     zion=nuclear charge 
c	(CAUTION:  I have only tested my transformation program for zion=1,
c 		   and I suspect it must be modified for other cases.      )
c	     (f,g)=(regular,irregular) energy normalized Coulomb functions.
c	     (fp,gp)= (d/dr)(f,g).
c	     Check:  W(f,g)=f*gp-fp*g  should be 2/pi if the program works.
c
c *********
c
      subroutine dummy
      implicit real*8(a-h,o-z)
      data pi/3.14159265358979d0/
      l=0
      eryd=-1.d0/(1.5d0**2) + 1.d-10
      zion=1.d0
      write(6,*) 'e=',eryd
      do 100 ir=1,10
        r=ir*0.5d0
        call seaton(l,eryd,r,zion,f,fp,g,gp)
        wfgtst=f*gp-fp*g - 2.d0/pi
        write(6,*) r,g,wfgtst
  100 continue
      epos=-eryd
      write(6,*) 'e=',epos
      do 200 ir=1,10
        r=ir*0.5d0
        call seaton(l,epos,r,zion,f,fp,g,gp)
        wfgtst=f*gp-fp*g - 2.d0/pi
        write(6,*) r,g,wfgtst
  200 continue
      stop
      end subroutine dummy
c
      subroutine seaton1(l,eryd,r,zion,f,fp,g,gp)
      implicit real*8(a-h,o-z)
      acc=1.d-10
      rl=l
      rr=r*zion
      eps=eryd/(zion**2)
      call coulfg(l,eps,rr,acc,f0,f0p,g0,g0p,k,ierr,actacc)
      if(.not.(eryd.lt.0))goto 23000
         ea=dabs(eryd)
         call ganda(a,gg,l,ea,zion,999)
         goto 23001
c     else
23000    continue
         gam=1.d0/dsqrt(eps)
         call gcmplx(a,gg,rl,gam)
23001 continue
      a5=dsqrt(dabs(a))
      f=a5*f0
      fp=a5*f0p
      g=(g0+gg*f0)/a5
      gp=(g0p+gg*f0p)/a5
c
c ** the next five lines changed on 1-22-88 by c.greene thanks to h. gao
c
	factor = dsqrt(zion)
	f=f/factor
	g=g/factor
	fp=fp*factor
	gp=gp*factor
c
      return
      end subroutine seaton1
c
c
      SUBROUTINE GANDA(A,G,L,E,ZION,NOPT)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA PI/3.1415926535897932D0/
      DPI = PI*2.D0
C*** THIS PROGRAM RETURNS THE QDT TRANSFORMATION PARAMETERS A&G.
      IF(ZION.EQ.0) WRITE(6,90)
   90 FORMAT(1X,'***** ZION = 0')
      IF(ZION.EQ.0) RETURN
      E = DABS(E)
      XNUE = ZION/DSQRT(E)
C*** EVALUATE A(K,L) FIRST.
      A = 1.D0
      IF(L.EQ.0) GO TO 109
      DO 100 I = 1,L
      A = A* (1.D0 -I*I*E/(ZION**2))
  100 CONTINUE
  109 continue
C*** GIVE WARNINGS IN CASE A < OR = 0 .
      IF(A.LE.0) WRITE(6,5555)
 5555 FORMAT(1X,'****** A < OR = 0')
      IF(NOPT.EQ.1) RETURN
C*** CHECK WHETHER XNUE = INTEGER.
      N = XNUE + 1.D-02
      Z = XNUE - N
      IF(Z.EQ.0) G = 0.D0
      IF(Z.EQ.0) RETURN
C*** G(K,L) IS NOW EVALUATED USING THE DIGAMMA PROGRAM.
      G = A*(DIGAM(L+1+XNUE) + DIGAM(-L+XNUE)
     1  - 2.D0 * DLOG(XNUE) )/DPI
      RETURN
      END subroutine ganda
      SUBROUTINE GCMPLX(B,G,RL,GAM)
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)
CCCC  COMPLEX*16 ZQ1,ZQ2,ZQ3,ZSUM,ZG
      COMPLEX*16 CDLOG
      DIMENSION GG(2)
      EQUIVALENCE (GG(1),ZG)
      DATA PI/3.1415926535897932D0/
C
      ZQ1 = DCMPLX(RL+1.D0,GAM)
      ZQ2 = DCMPLX(-RL,GAM)
      ZQ3 = DCMPLX(0.D0,GAM)
      ZSUM = ZDIGAM(ZQ1) + ZDIGAM(ZQ2) + CDLOG(ZQ3)*(-2.D0)
      PROD = 1.D0
      L = RL
      DO 100 I = 1,L
      IF(L.EQ.0) GO TO 100
      PROD = PROD*( 1.D0 + I*I/(GAM*GAM) )
  100 CONTINUE
      A1 = PROD
      ZG = ZSUM*A1/(2.D0*PI)
      G = GG(1)
      QQ = 1.D0 - DEXP(-2.D0*PI*GAM)
      B = A1 / QQ
      RETURN
      END subroutine gcmplx
c
c      function cdabs(z)
c      implicit real*8(a-h,o-y),complex*16(z)
c      cdabs = cabs(z)
c      return
c      end
ccc
c	function cdlog(z)
c	implicit complex*16(a-h,o-z)
c	cdlog=clog(z)
c	return
c	end
c
c      function cdexp(z)
c      implicit complex*16(a-h,o-z)
c      cdexp=cexp(z)
c      return
c      end
c
      FUNCTION DIGAM(ARGG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIGAM = 0.D0
      ARG5 =-0.57721566490153286D0
      EN = 1.D0
      ARG  = ARGG
      ARG2 = ARGG
    1 IF(ARG2-40.D0) 2,3,3
    2 DIGAM = DIGAM - 1.D0/ARG
      ARG = 1.D0+ ARG
      ARG2 = 1.D0 + ARG2
      GO TO 1
    3 PSI = DLOG(ARG) - 1.D0/(2.D0*ARG) - 1.D0/(12.D0*ARG**2)
      PSI = PSI + 1.D0/(120.D0*ARG**4) - 1.D0/(252.D0*ARG**6)
      DIGAM = DIGAM + PSI
      RETURN
      END function digam

      FUNCTION ZDIGAM(arg)
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)
      COMPLEX*16 ARG,ZDIGAM,cdlog
      COMPLEX*8 ARGG
      REAL*8 INC
      ZDIGAM = (0.D0,0.D0)
      ARG5 =-0.57721566490153286D0
      PI = 3.1415926535897932D0
      EN = 1.D0
      ARGG = ARG
      ARG2 = REAL(ARGG)
      ARG3 = AIMAG(ARGG)
      IF(ARG3) 4,1,4
    1 IF(ARG2-40.D0) 2,3,3
    2 ZDIGAM = ZDIGAM - 1.D0/ARG
      ARG = 1.D0+ ARG
      ARG2 = 1.D0 + ARG2
      GO TO 1
    3 PSI = CDLOG(ARG)-1.D0/(2.D0*ARG)-1.D0/(12.D0*ARG**2)
      PSI=PSI +1.D0/(120.D0*ARG**4)-1.D0/(252.D0*ARG**6)
      ZDIGAM = ZDIGAM + PSI
      GO TO 12
    4 IF(ARG2) 5,7,6
    5 ZDIGAM = ZDIGAM - 1.D0/ARG
      ARG = ARG + 1.D0
      ARG2 = ARG2 + 1.D0
      GO TO 4
    6 ARG = ARG - 1.D0
      ARG2 = ARG2 - 1.D0
      ZDIGAM = ZDIGAM + 1.D0/ARG
      GO TO 4
    7 Y = CDABS(ARG)
      ARG7 = PI*Y
      ARG4 = 0.5D0/Y + (PI/2.D0)/DTANH(ARG7)
      IF(Y-20.D0) 8,10,10
    8 INC = Y*Y/(EN*(EN*EN+Y*Y))
      ARG5 = ARG5 + INC
      IF(INC - 1.D-12) 11,11,9
    9 EN = EN + 1.D0
      GO TO 8
   10 ARG5 = 1.D0/(12.D0*Y**2) + 1.D0/(120.D0*Y**4)
      ARG5 = ARG5 + 1.D0/(252.D0*Y**6)+ DLOG(Y)
   11 ZDIGAM = DCMPLX(ARG5,ARG4) + ZDIGAM
C     XQ1 = REAL(ZDIGAM)
C     XQ2 = AIMAG(ZDIGAM)
   12 continue
      return
      END function zdigam
	subroutine seaton(l,eryd,r,zion,f,fp,g,gp)
	implicit real*8(a-h,o-z)
c -- assume for now that zion=1
	rho=zion*r
	epsr=eryd/(zion**2)
	acc=1.d-11
	call fogo(f,fp,g,gp,epsr,l,rho,ww,acc,actacc)
	s2=dsqrt(2.d0)
	f=f*s2
	fp=fp*s2
	g=g*s2
	gp=gp*s2
chg -- add the following lines on 1-8-91:

	factor = dsqrt(zion)
	f=f/factor
	g=g/factor
	fp=fp*factor
	gp=gp*factor

chg -- end additions

	return
	end subroutine seaton
cccc
      SUBROUTINE FOGO(F0,FP0,G0,GP0,EPSR,LL,RO,WW,ACC,ACTACC)
*
*   INPUT EPSR (E RYDBERG) LL (L) ACC=1.E-11
*    OUTPUT ENERGY NORMALISED FUNCTIONS AND DERIVATIVES (RYDBERG)
*
*
        IMPLICIT REAL*8 (A-H,O-Z)
chg
	common/agqdt/aqdt,gqdt,fr0,fpr0,gr0,gpr0
chg
      DATA R2PI/.159154943 09189 5336D0/
       DATA PI/ 3.1415 92653 58979 32D0/
C     PRINT 1013
1013  FORMAT(//////,121(1H*) ,/)
*
       A=1.D0
       IF(LL.EQ.0) GO TO 50
       TESEN=1.D-10
       IF(DABS(EPSR).LT.TESEN)  GO TO 50
C
C
C CALCUL DE A= PRODUIT( 1. + EPS * P2 )
       DO 20 LP=1,LL
       B=1.D0 + LP*LP*EPSR
       A=A*B
C     PRINT *,' A B ',A,B
20     CONTINUE
50     CONTINUE
C     PRINT 1010 ,EPSR,LL,RO,A
1010  FORMAT(/////,' ENERGIE  ',D25.16,' RYD  L= ',I3,' RHO = ',
     =D25.16,/,' A = ',D25.16 ,//)
C
C
C
C CALCUL DE G ET DE B
C *******************
      IF(DABS(EPSR).LT.TESEN )  GO TO 8000
      IF(EPSR.GE.TESEN) GO TO 5000
C
C CAS ENERGIE NEGATIVE
      PNU=DSQRT(-1.D0/EPSR)
      PL1=DFLOAT(LL)
chg
	pl1=0.d0
chg
      IF(PNU.GT.PL1)  GO TO 61
      PRINT 1062,LL,EPSR
 1062 FORMAT(//////,131(1H*),/,
     1'****** CAS ETUDIE IMPOSSIBLE  .  L = ',
     2I4,' ENERGIE =  ',D25.13,' RYDBERG ',/, 131(1H*) )
      STOP
 61   CONTINUE
      X=PNU+LL+1.D0
      PSI=D1IGAM(X)
C     PRINT *,' X PSI ',X,PSI
      X=PNU-LL
      PSI=PSI +D1IGAM(X) -2.D0*DLOG(PNU)
      X2=DLOG(PNU)
      X1=D1IGAM(X)
C     PRINT *,' X PSI LOG ',X,X1,X2
      G=A*R2PI*PSI
C     PRINT 1020,G,PNU,PSI
1020  FORMAT('  G = ',D25.16,' KAPA ',D25.16,' PSI',D25.16)
      B=A
C     PRINT *, ' B = ',B
chg
	aqdt=b
	gqdt=g
chg
      GO TO 4000
5000    CONTINUE
C
C CAS ENERGIE POSITIVE
      PGA=DSQRT(1.D0/EPSR)
      SUM=0.D0
C CALCUL DE LA DEPENDANCE EN L
      IF(LL.EQ.0) GO TO 90
      DO 80 LP=1,LL
      S1=LP*EPSR
      SUM=SUM + S1/(1.D0 + LP*S1)
C     PRINT *,'  SUM S1  ',SUM,S1
 80   CONTINUE
90    CONTINUE
      IF(EPSR . GT. 0.05D0 ) GO TO 92
C ENERGIE < 0.05 RYD . FORMULE 2
      X1=GAMI(PGA,2)
      GO TO 99
 92   IF(EPSR . GT. 0.8D0 ) GO TO 96
C 0.05 RYD < ENERGIE < 0.8 RYD . FORMULE 1
      X1=GAMI(PGA,1)
      GO TO 99
C O.8 RYD < ENERGIE . FORMULE 3
96    X1=GAMI(PGA,3)
99    CONTINUE
      G=(SUM + X1)* A * 2.D0 * R2PI
C     PRINT 1040,G,PGA,SUM,X1
 1040 FORMAT('  G = ',D25.16,' GAMA = ',D25.16,' SUM = ',
     +D25.16,' X = ',D25.16)
      B=A/(1.D0 - DEXP( -PGA/R2PI))
C     PRINT *, '  B = ',B
chg
	aqdt=b
	gqdt=g
chg
      GO TO 4000
8000   CONTINUE
C
C CAS ENERGIE=0
      G=0.D0
      B=1.D0
C     PRINT 1050,G
1050  FORMAT('  G = ',D25.16)
C     PRINT *,'  B = ',B
4000  CONTINUE
C
C CALCUL DE F(R0) ET DE G(R0) ANALYTIQUES EN ENERGIE
C **************************************************
C*********WRONSKIEN = 2 / PI
C         K TERMES DANS LE CALCUL
C         ACTACC PRECISION RELATIVE DU WRONSKIEN
C
      CALL COULFG(LL,EPSR,RO,ACC,FR0,FPR0,GR0,GPR0,KK
     1,IERR,ACTACC)
C     PRINT *
C     PRINT *
C     PRINT *,'  LL  EPSR   RO ',LL,EPSR,RO
C     PRINT *,' VALEURS DE F  FP  G  GP '
C     PRINT *,FR0,FPR0,GR0,GPR0
C     PRINT *,' K  IER WRONSKIEN ',KK,IERR,ACTACC
C
C CALCUL DE F0 ET G0 NORMEES EN RYDBERG
C *************************************
C
chg
	b=dabs(b)
chg
      HR0 =GR0 + G * FR0
      HPR0 = GPR0 + G * FPR0
      B = DSQRT ( B / 2.D0  )
      F0 = FR0 * B
      FP0 = FPR0 * B
      G0 = HR0 / ( 2.D0 * B )
       GP0 = HPR0 /( 2.D0 * B )
      WW=(F0*GP0-FP0*G0)*PI-1.D0
C       PRINT*,' FONCTION NORMEES EN RYDBERG'
C     PRINT *,F0,FP0,G0,GP0
C     PRINT*,' 2EME WRONSKIEN ',WW
C      PRINT *
      RETURN
      END subroutine fogo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION GAMI(X,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(7)
      DIMENSION C(20)
      DIMENSION CCOR(4)
      DATA CCOR/ .20200 74006 59677 52D0,.03692 77526 92953 20D0,
     1 .0083492773817610 2D0, .00200839282608212D0 /
      DATA B /8.33333 33333 33333D-2,+8.33333 33333 33333D-3,
     1        3.96825 39682 53968D-3,+4.16666 66666 66667D-3,
     2        7.57575 7575757 576D-3,+2.10927 96092 79609D-2,
     3        8.33333 33333 33333D-2/
      DATA C /2.02056 90315 95943D-1, 3.69277 55143 36993D-2,
     1       8.34927 73819 22827D-3, 2.00839 28260 82214D-3,
     2       4.94188 60411 94646D-4, 1.22713 34757 84891D-4,
     3       3.05882 36307 02049D-5, 7.63719 76378 9976 D-6,
     4     1.90821 27165 5394 D-6, 4.76932 98678 781  D-7,
     5       1.19219 92596 531  D-7, 2.98035 0351465  D-8,
     6       7.45071 17898 4   D-9, 1.86265 97235 1 D-9,
     7       4.65662 90650    D-10, 1.16415 50173  D-10,
     8  2.910385044 D-11, 7.27595 984 D-12,
     9       1.81898 965    D-12,4.547473 8 D-13 /
      DATA EUL/.5772 15664 90153  29 D0/
      DATA TES/1.D-10/
C CALCUL DE REEL(PSI(I*GAM)) - LOG(GAM)
C        N=1 SOMME 1/N(N2+GAM2)
C        N=2 DEVELOPPEMENT ASYMPTOTIQUE
C        N=3 CAS GAM INF A 2
      IF(N.GE.4) GO TO 8000
      A=DABS(X)
      AA=A*A
      IF(N.EQ.2) GO TO 50
      IF(N.EQ.3) GO TO 100
      S1=1.D0/(1.D0 + AA)
C     PRINT *, '  S1  ',S1
      DO 10 IK=2,100
      DS=AA + IK*IK
      DS=1.D0/(DS*IK)
      S1=S1+DS
10    CONTINUE
C     PRINT *,' X S1 DS TES ',X,S1,DS,TES
C     PRINT *,' ***** DEVELOPPEMENT NON CONVERGE'
      GO TO 30
20    S1=S1+DS
30    CONTINUE
      RAP=DS/S1
      S2=DLOG(A)
C     PRINT *,' S1,LOG',S1,S2
      GAMI=S1*AA - EUL -S2
       RAP2=DABS(GAMI/(EUL+S2))
C     PRINT *,' CONVERG ANULATION ',RAP,RAP2
C     PRINT 1020,N,X,GAMI
1020  FORMAT(' OPTION N= ',I3,' X = ',D25.16,' GAMI = ',
     1D25.16,/)
      C1=0.D0
      X1=-X*X
      X2=-1.D0
      DO 5010 IC=1,4
      X2=X2*X1
      C2=X2*(C(IC)  - CCOR(IC))
      C1=C1 + C2
      GAMI= GAMI + C2
C     PRINT *,' C1 C2 GAMI ',C1,C2,GAMI
5010  CONTINUE
C     PRINT 1020,N,X,GAMI
      RETURN
50    AA=1.D0/AA
      A1=1.D0
      S1=0.D0
      DO 60 IK=1,7
      A1=A1*AA
      DS=B(IK)*A1
      S1= S1 + DS
      RAP=DS/S1
C     PRINT *,' S1 DS RAP ',S1,DS,RAP
60    CONTINUE
      GAMI=S1
C     PRINT 1020,N,X,GAMI
      RETURN
100   A1=-1.D0
      S1=0.D0
      DO 110 IK=1,20
      A1=-A1*AA
      DS=C(IK)*A1
      S1=S1 + DS
      RAP=DS/S1
C     PRINT *,' S1 DS RAP ', S1,DS,RAP
110   CONTINUE
      S2=DLOG(A)
      S3=1.D0/(1.D0 + AA)
      GAMI = 1.D0 -EUL -S2+S1 -S3
C     PRINT *,' S1 LOG S3 ',S1,S2,S3
      RAP=DABS(S1/GAMI)
      RAP=1.D0/RAP
C     PRINT *,'  ANULATION ',RAP
C     PRINT 1020,N,X,GAMI
      RETURN
8000   PRINT 8010,N
      GAMI=0.D0
8010      FORMAT('  FORMULE NON PROGRAMMEE N',
     +' ***********************',///)
      RETURN
       END function gami
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION D1IGAM(X)
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION B(6)
      DATA PI/ 3.14159 26535 89793 D0/
      DATA B / +8.33333 33333 33333D-2, -8.33333 33333 33333D-3,
     1         +3.96825 39682 53968D-3, -4.16666 66666 66667D-3,
     2         +7.57575 75757 57576D-3, -2.10927 96092 79609D-2/
C CALCUL DE PSI(X)  X REEL
C PASSAGE X POSITIF PLUS GRAND QUE 15
C FORMULE ASYMPTOTIQUE A 7 TERMES
      A=DABS(X)
C     IF(DINT(X) + A)  ,4,
      XZ=DINT(X)+A
      IF(XZ.EQ.0.D0) GO TO 4
      V=A
      H=0.D0
      IF(A.GE.(15.D0)) GO TO 3
      N=14 - INT(A)
      H=1.D0/V
C     IF(N)  ,2,
      IF(N.EQ.0) GO TO 2
      DO 1 I=1,N
      V= V + 1.D0
1     H= H + 1.D0/V
2     V= V + 1.D0
3     R=1.D0/V**2
      D1IGAM=DLOG(V) -0.5D0/V -R*(B(1)+R*(B(2)+R*(B(3)+R*
     1 (B(4)+R*(B(5)+R*(B(6)+R*B(1))))))) - H
      IF(X . GE.(0.000D0) ) RETURN
      H=PI*A
      D1IGAM=D1IGAM + 1.D0/A + PI*DCOS(H)/DSIN(H)
      RETURN
C SORTIE ERREUR  :  X ENTIER NON POSITIF
4     PRINT 100,X
100   FORMAT(/////,131(1H*)//, '  *** DDIGAM
     1  ARGUMENT ENTIER NON NEGATIF =',D16.9,
     2  ' ***** ' )
      D1IGAM=0.D0
      RETURN
      END function d1igam

      end module gensub
