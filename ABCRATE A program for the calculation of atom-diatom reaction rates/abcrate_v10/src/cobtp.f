!!!*************************************************************
! 文件/File: cobtp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cobtp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE COBTP (E,Q1,Q2,NOM,FB,AB,GB)
C
C     COBTP  - find turning points in centrifugal oscillator potential
C
C  Called by:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     COBINT - compute phase integrals needed for WKB quantization
C
C  Calls:
C     CUBIC  -  computes roots of cubic equation   
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION RRT(3),AIRT(3)
C
C      write (6,*) ' cobtp called with nom,fb,ab,gb,e='
C      write (6,*) nom,fb,ab,gb,e
      IF (NOM.EQ.0) THEN
C  omega=0 roots
         IF (AB.EQ.0.0D0) THEN                                          GCL1092
C            write (6,*) ' ab=0'
            Q1 = 0.0D0                                                  GCL1092
            Q2 = SQRT(2.0D0*E/FB)                                       GCL1092
         ELSE IF (FB.EQ.0.0D0) THEN                                     GCL1092
C            write (6,*) ' fb=0'
            Q1 = 0.0D0                                                  GCL1092
            Q2 = (24.0D0*E/AB)**0.25D0                                  GCL1092
         ELSE 
C  solve quadratic equation for roots
            T1 = -6.0D0*FB/AB                                           GCL1092
            T2 = SQRT(1.0D0 + E*AB/(1.5D0*FB*FB))                       GCL1092
C            write (6,*) ' t1,t2=', t1,t2
            IF (FB.LT.0.0D0) THEN                                       GCL1092
               Q2 = SQRT(T1*(1.0D0+T2))                                 GCL1092
C               write (6,*) ' fb<0, q2=', q2
               IF (E.LT.0.D0) THEN                                      GCL1092
                  Q1 = SQRT(T1*(1.0D0-T2))                              GCL1092
               ELSE
                  Q1 = 0.0D0                                            GCL1092
               END IF
            ELSE
               Q1 = 0.0D0                                               GCL1092
               Q2 = SQRT(T1*(1.0D0-T2))                                 GCL1092
            END IF
         END IF
C
      ELSE
C  omega#0, need to solve cubic equation for roots
         A = AB/24.0D0                                                  GCL1092
         B = 0.5D0*FB                                                   GCL1092
         C = -E
         D = 0.5D0*DBLE(NOM*NOM)*GB                                     GCL1092
c         write (6,*) ' in cobtp, a,b,c,d='
c         write (6,*) a,b,c,d
         CALL CUBIC(A,B,C,D,NREAL,RRT,AIRT)
c         write (6,*) ' nreal,roots=', nreal
c         write (6,*) (rrt(i),airt(i),i=1,3)
         NRT = 0
C  check that real roots are left or right bounds
         IF (NREAL.GT.0) THEN
            AD = AB/6.0D0                                               GCL1092
            BD = FB
            CD = -2.0D0*D                                               GCL1092
            Q1 = 1000.0D0                                               GCL1092
            Q2 = -1000.0D0                                              GCL1092
            DO 10 I = 1,NREAL
               Q0SQ = RRT(I)
               NRT = NRT + 1
               IF (Q0SQ.GT.0.0D0) THEN                                  GCL1092
                  Q0 = SQRT(Q0SQ)
                  DV = Q0*(AD*Q0SQ + BD) + CD/(Q0*Q0SQ)
                  IF (DV.LT.0.0D0) THEN                                 GCL1092
                     Q1 = Q0
                  ELSE
                     Q2 = Q0
                  END IF
               END IF
10          CONTINUE
         END IF
         IF (NRT.LT.2 .OR. Q1.GT.Q2) THEN
            WRITE (6, 6000)
            STOP 'COBTP 1'
         END IF
      END IF
      RETURN
6000  FORMAT(/,2X,T5,'Error: Problem finding the turning points ',
     *               'for the centrifugal oscillator.')
      END
