!!!*************************************************************
! 文件/File: cobvex.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cobvex.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE COBVEX (NOM,FB,AB,GB,V0,VMAX)
C
C     COBVEX   - find extrema in centrifugal oscillator potential
C
C  Call by:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C
C  Calls:
C     CUBIC  -  computes roots of cubic equation   
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION RRT(3),AIRT(3)
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
C
      LMAX = .FALSE.
      VMAX = 1.0D6                                                      GCL1092
      IF (FB.LE.0.0D0 .AND. AB.LE.0.0D0) THEN                           GCL1092
         LMAX = .TRUE.
         RETURN
      ELSE IF (NOM.EQ.0) THEN
C  omega = 0
         IF (FB.GT.0.0D0) THEN                                          GCL1092
            V0 = 0.0D0                                                  GCL1092
            IF (AB.LT.0.0D0) THEN                                       GCL1092
C    for AB<0, maximum in potential
               Q0SQ = -6.0D0*FB/AB                                      GCL1092
               VMAX = -1.5D0*FB*FB/AB                                   GCL1092
            END IF
         ELSE
C    for FB<0,AB>0, minimum in potential
            Q0SQ = -6.0D0*FB/AB                                         GCL1092
            V0= -1.5D0*FB*FB/AB                                         GCL1092
         END IF
      ELSE IF (AB.EQ.0.0D0) THEN                                        GCL1092
C  omega#0, AB=0, minimum 
         Q0SQ = SQRT(GB/FB)*DBLE(NOM)                                   TCA0497
         V0 = FB*Q0SQ
      ELSE IF (FB.EQ.0.0D0) THEN                                        GCL1092
C  omega#0, FB=0, AB>0, minimum
         Q0SQ = (6.0D0*GB*DBLE(NOM*NOM)/AB)**(1.D0/3.D0)                TCA1196
         V0 = AB*Q0SQ*Q0SQ/8.0D0                                        GCL1092
      ELSE
C  omega#0, solve cubic equation for roots (squares of extrema)
         A = AB/6.0D0                                                   GCL1092
         B = FB
         C = 0.0D0                                                      GCL1092
         D = -DBLE(NOM*NOM)*GB                                          TCA0497
C         write (6,*) ' in cobvex, a,b,c,d='
C         write (6,*) a,b,c,d
         CALL CUBIC(A,B,C,D,NREAL,RRT,AIRT)
C         write (6,*) ' nreal=', nreal
C         do 6666 i = 1,3
C            t = d + rrt(i)*(c + rrt(i)*(b + rrt(i)*a))
C            write (6,*) ' root', i, 'f(root)=', rrt(i),airt(i),t
C6666     continue
         NMN = 0
         NMX = 0
         IF (NREAL.GT.0) THEN
C  for real roots, check if minimum or maximum
            A = AB/24.0D0                                               GCL1092
            B = 0.5D0*FB                                                GCL1092
            C = 0.5D0*DBLE(NOM*NOM)*GB                                  GCL1092
            AD = 0.5D0*AB                                               GCL1092
            BD = FB
            CD = 6.0D0*C                                                GCL1092
            DO 10 I = 1,NREAL
               Q0SQ = RRT(I)
               IF (Q0SQ.GT.0.0D0) THEN                                  GCL1092
                  D2V = AD*Q0SQ + BD + CD/(Q0SQ*Q0SQ)
                  V = Q0SQ*(A*Q0SQ + B) + C/Q0SQ
C                  write (6,*) ' i,q0sq,v,d2v='
C                  write (6,*) i,q0sq,v,d2v
                  IF (D2V.LT.0.0D0) THEN                                GCL1092
                     NMX = NMX + 1
                     IF (NMX.GT.1) THEN
                        VMAX = MAX(V,VMAX)
                     ELSE
                        VMAX = V
                     END IF
                  ELSE
                     NMN = NMN + 1
                     IF (NMN.GT.1) THEN
                        V0 = MIN(V0,V)
                     ELSE
                        V0 = V
                     END IF
                  END IF
               END IF
10          CONTINUE
         END IF
         IF (NMN.EQ.0 .AND. NMX.EQ.0) LMAX = .TRUE.
      END IF
      RETURN
      END
