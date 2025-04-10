!!!*************************************************************
! 文件/File: call_bes.sub.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: call_bes.sub.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine call_bessel_complex(x,lambda,charge,rad,fc_ret,gc_ret,fcp_ret,gcp_ret,kfn)
  !  program call_coulomb_complex
  use nrtype, only : dbl,dpc,i4b
  implicit none
  real(kind=dbl) :: epsi,zero,half,one,four,Pi,rerr
  real(kind=dbl),intent(in) :: charge,rad
  complex(kind=dpc),intent(in) :: lambda
  COMPLEX(kind=dpc) :: X,ETA,ZLMIN,FC(201),GC(201),FCP(201),GCP(201), &       
       & SIG(201),ZL,WS,CI,CGAMMA
  complex(kind=dpc),intent(out) :: fc_ret,gc_ret,fcp_ret,gcp_ret                        
  INTEGER(kind=i4b),intent(in) :: kfn 
  INTEGER(kind=i4b) LDISP(8),nl,mode,ifail,md,ih,kfin,ii,jj, &
       & nfp,npq,kase,i,l,ipr                                                 
  LOGICAL WHIT                                                     
  COMMON       /STEED/ RERR,NFP(2),NPQ(3),KASE(2)                   
  CHARACTER*20 NOTE                                                 
  CHARACTER*4 WHO(2,4,3),IRREG,REG                                  
  DATA ZERO,HALF,ONE,FOUR,CI / 0D+0,0.5D+0,1D+0,4D+0,(0D+0,1D+0) /  
  DATA WHO / 'F','G','j','y','J','Y','I','K' ,           &           
       &           'F','H+','j','h(1)','J','H(1)','?','?' ,    &           
       &           'F','H-','j','h(2)','J','H(2)','?','?' /               
  DATA LDISP / 1,2,4,11,31,101,301,1001 /                           
  !                                                                       
  ipr=2 ! Print error messages
  PI = FOUR * ATAN(ONE)                                            
  !WRITE(6,1000)                                                     
  !                                                                       
10 continue

  !-----------------------------------------------------------------------
  ! Define input
  !-----
  note='mostreo'
  !READ(5,*,END=40) X,ETA,ZLMIN,NL,MODE,KFN,WHIT,NOTE
  !ene=.002d0
  !charge=1.d0
  !rad=8.d0
  !lambda=(0.d0,0.d0)
  !-----
  mode=1
  nl=1
  whit=.FALSE.

  !-----------------------------------------------
  ! Make the variables complex
  !  x=cmplx(Sqrt(2.d0*Ene)*rad,0.d0)
  !  eta=cmplx(-charge/Sqrt(2.d0*Ene),0.d0)
  !x= sqrt(cmplx(2.d0*Ene))*rad
  !eta=-charge/sqrt(cmplx(2.d0*Ene))
  !------------------------------------------------

  zlmin=lambda
  !write(3004,*)x,eta,zlmin,nl,mode,kfn,ene,charge,lambda,rad,RERR,NFP,NPQ,KASE
  !-----------------------------------------------------------------------
  IF(NL.LE.0) GO TO 40                                              
  IFAIL = 1                                                         
  MD = MOD(ABS(MODE),10)                                            
  IH =     ABS(MODE)/10                                             
  KFIN = KFN                                                        
  !WRITE(6,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE                        
  !IF(WHIT) KFN=MAX(KFN,0)                                           
  !WRITE(6,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE                        
  !                                                                       
  !do jj =1, 3000 
  !   epsi=4.d0
  !   eta=cmplx(-1.d0/sqrt(epsi),0.d0)
  !   x=(0.1d0+(jj-1)*.05d0)  * 2.d0
  !   ii=0
  CALL WCLBES(X,ETA,ZLMIN,NL,FC,GC,FCP,GCP,SIG,KFN,MODE,IFAIL,IPR)      
  !write(6,*)X,ETA,ZLMIN,NL,FC(1),GC(1),FCP(1),GCP(1),SIG(1),MODE,KFN,IFAIL
  !write(3004,*)x,eta,zlmin,nl,mode,kfn,ene,charge,lambda,rad,RERR,NFP,NPQ,KASE
  !     Conversion into Seaton's form
  do ii=0,nl-1
  !   fc(ii+1)=fc(ii+1)*((-2.d0*eta)**(zlmin+ii+1)*exp(Pi*eta/2.d0) /(2.d0**  &
  !        &  (zlmin+ii)*abs(cgamma(zlmin+ii+1.d0+eta*cmplx(0.d0,1.d0)))))
  !   gc(ii+1)=gc(ii+1)*2.d0**(zlmin+ii)*abs(cgamma(zlmin+ii+1.d0+eta*cmplx(0.d0,1.d0)))&
  !        &*exp(-Pi* eta/2.d0)/(Pi*(-2.d0*eta)**(zlmin+ii)) 
  !   fcp(ii+1)=fcp(ii+1)*((-2.d0*eta)**(zlmin+ii+1)*exp(Pi*eta/2.d0) /(2.d0** &
  !        &  (zlmin+ii)*abs(cgamma(zlmin+ii+1.d0+eta*cmplx(0.d0,1.d0))))) / (charge/sqrt(2.d0*ene))
  !   gcp(ii+1)=gcp(ii+1)*2.d0**(zlmin+ii)*abs(cgamma(zlmin+ii+1.d0+eta*cmplx(0.d0,1.d0)))&
  !        &*exp(-Pi* eta/2.d0)/(Pi*(-2.d0*eta)**(zlmin+ii)) / (charge/sqrt(2.d0*ene))

     !      fc(ii+1)=fc(ii+1)*(2*exp((eta*Pi)/2.)*(-eta)**(1 + l))
     !     & /Abs(CGamma(1.d0 +(0,1)*eta + l))
  end do
  !write(35,*)real(x)/2.d0,gcp(1)
  !write(45,*)real(x)/2.d0,fc(1)
  !write(55,*)real(x)/2.d0,fcp(1)
  !write(65,*)real(x)/2.d0,gc(1)
  !  ,cgamma(1+eta*(0,1))
  !write(34,*)real(x)/2.d0,real(gc(2))
  !write(33,*)real(x)/2.d0,real(gc(3))
  !write(32,*)real(x)/2.d0,real(gc(4))
  !end do
  !                                                                       
  !do ii=1,nl
  !     Conversiion into Seaton's form
  !      fc(ii)=fc(ii)*(-2.d0*eta)**(zlmin+ii+1)*exp(Pi*eta/2.d0)/(2.d0** 
  !     &  (zlmin+ii)*cgamma(zlmin+ii+1+eta*cmplx(0.d0,1.d0)))
  !      write(35,*)real(x),real(fc(1))
  !                                                                      
  !write(36,*)fc(ii),gc(ii),fcp(ii),gcp(ii),x,eta,zlmin,nl,ii
  !end do
  fc_ret=(fc(1))
  gc_ret=(gc(1))
  fcp_ret=(fcp(1))
  gcp_ret=(gcp(1))
  !write(307,*)fc_ret,gc_ret,fcp_ret,gcp_ret
  !WRITE(6,1020) IFAIL,RERR,NFP,NPQ,KASE                             
  !IF(IFAIL.LT.0) GO TO 30                                           
  !DO 20 I=1,8                                                       
  !   L = LDISP(I)                                                     
  !   IF(L.GT.NL-IFAIL) GO TO 20                                        
  !   ZL = ZLMIN + L - 1                                             
  !   IF(KFN.NE.0) SIG(L) = ZERO                                     
  !   IRREG = WHO(2,MAX(KFN+1,1),IH+1)                               
  !   REG = WHO(1,MAX(KFN+1,1),1)                                  
  !   IF(WHIT) THEN                                                  
  !      IRREG = 'WHIT'                                              
  !      WS = EXP(-HALF*PI*(ETA - CI*ZL)  - CI*SIG(L))               
  !      GC(L)  = WS * GC(L)                                         
  !      IF(MD.EQ.1) GCP(L) = CI*WS * GCP(L)                         
  !      FC(L)  = CI/WS * FC(L)                                      
  !      IF(MOD(MD,2).EQ.1) FCP(L)  = FCP(L) / WS                    
  !      REG = 'WH-F'                                               
  !   ENDIF
  !   !WRITE(6,1030) ZL,REG,FC(L),IRREG,GC(L)                            
  !   IF(MD.EQ.1)WRITE(6,1040)            FCP(L),GCP(L)                 
  !   IF(SIG(L).NE.ZERO.AND.KFIN.EQ.0) WRITE(6,1050) SIG(L)             
20 CONTINUE                                                          
30 CONTINUE                                                          
  !     GO TO 10                                                          
40 continue
  !STOP                                                              
1000 FORMAT('1TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL ROUTINES'&
       & /)                                                               
1010 FORMAT(/'0X =',F12.4,F9.4,', ETA =',2F8.3,', ZLMIN =',2F8.3,'  NL &
       &=',I4,  '  MODE = ',I3,'  KFN =',I2,8X,A20)                       
1020 FORMAT(' COULCC  :: IFAIL =',I4,                                  &
       &   '  RERR =',1P,E10.1,'  ITS =',I6,I4,2I6,I4,2I2/)               
1030 FORMAT(' ZL =',2F8.3,' :: FC = ',A4,' =',1P,2D20.12,',  GC = ',   &
       &   A4,' =',2D20.12)                                               
1040 FORMAT(24X,' FC''=       ',1P,2D20.12,',  GC'' =',6X ,2D20.12)    
1050 FORMAT(24X,' SIG=   ',   2F20.12)                                 
end subroutine call_bessel_complex
!end program call_coulomb_complex
