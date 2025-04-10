!!!*************************************************************
! 文件/File: ellip.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: ellip.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE ELLIP (X, EL1, EL2)
C
C     ELLIP  - elliptic integrals of first and second kind
C
C  Called by:
C     HQSC   - compute semiclassical eigenvalues of harm.-quart.
C     THETA2 - phase integral for semiclassical eigenvalues of double
C              well
C     THETA3 - phase integral for semiclassical eigenvalues of double
C              well
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(20), SUM(4)
      SAVE A                                                            TCA1097
      DATA A/
     *   1.451196212D-2,  4.41787012D-3, 1.736506451D-2,  5.26449639D-3,
     *   3.742563713D-2, 3.328355346D-2, 4.757383546D-2, 4.069697526D-2,
     *   3.590092383D-2, 6.880248576D-2, 6.260601220D-2, 9.200180037D-2,
     *   9.666344259D-2,1.2498593597D-1,4.4325141463D-1,2.4998368310D-1,
     * 1.38629436112D00,         0.5D00,         1.0D00,         0.0D00/
C
C  elliptic integrals of the first and second kind
      IF (X .LT. 0.0D0 .OR. X .GE. 1.0D0) THEN                          GCL1092
         WRITE (6, 6000)
         STOP 'ELLIP 1'
      END IF
      XP = 1.0D0 - X
      DO 10 I = 1,4
         SUM(I) = A(I)
   10 CONTINUE
      K = 4
      DO 30 J = 1,4
         DO 20 I = 1,4
            K = K + 1
            SUM(I) = A(K)+SUM(I)*XP
   20    CONTINUE
   30 CONTINUE
      T = -LOG(XP)
      EL1 = SUM(1) + T*SUM(2)
      EL2 = SUM(3) + T*SUM(4)
      RETURN
6000  FORMAT(/,2X,T5,'Error: In ELLIP X must be .GE. 0 .AND. .LT. 1')
      END
