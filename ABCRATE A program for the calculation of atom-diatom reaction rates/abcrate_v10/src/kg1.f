!!!*************************************************************
! 文件/File: kg1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: kg1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE KG1 (N, NDIM, PT, WT)
C
C     KG1    - set up Kronrod quadrature points and weights
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION PT(NDIM), WT(NDIM,2)
      M = 2*N + 1
      IF (N .EQ. 2) THEN
C         (N,2N+1)=    2    5
         PT( 1) = -9.25820099772551D-01
         PT( 2) = -5.77350269189619D-01
         PT( 3) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  1.00000000000000D+00
         WT( 3, 1) =                 0.D0
         WT( 1, 2) =  1.97979797979798D-01
         WT( 2, 2) =  4.90909090909090D-01
         WT( 3, 2) =  6.22222222222224D-01
      ELSE IF (N .EQ. 3) THEN
C         (N,2N+1)=    3    7
         PT( 1) = -9.60491268708019D-01
         PT( 2) = -7.74596669241468D-01
         PT( 3) = -4.34243749346798D-01
         PT( 4) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  5.55555555555539D-01
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  8.88888888888832D-01
         WT( 1, 2) =  1.04656226026467D-01
         WT( 2, 2) =  2.68488089868333D-01
         WT( 3, 2) =  4.01397414775962D-01
         WT( 4, 2) =  4.50916538658474D-01
      ELSE IF (N .EQ. 4) THEN
C         (N,2N+1)=    4    9
         PT( 1) = -9.76560250737570D-01
         PT( 2) = -8.61136311594013D-01
         PT( 3) = -6.40286217496310D-01
         PT( 4) = -3.39981043584844D-01
         PT( 5) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  3.47854845137421D-01
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  6.52145154862517D-01
         WT( 5, 1) =                 0.D0
         WT( 1, 2) =  6.29773736654728D-02
         WT( 2, 2) =  1.70053605335723D-01
         WT( 3, 2) =  2.66798340452285D-01
         WT( 4, 2) =  3.26949189601452D-01
         WT( 5, 2) =  3.46442981890137D-01
      ELSE IF (N .EQ. 5) THEN
C            (N,2N+1)=    5   11
         PT( 1) = -9.84085360094845D-01
         PT( 2) = -9.06179845938627D-01
         PT( 3) = -7.54166726570851D-01
         PT( 4) = -5.38469310105661D-01
         PT( 5) = -2.79630413161783D-01
         PT( 6) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  2.36926885056196D-01
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  4.78628670499345D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  5.68888888888814D-01
         WT( 1, 2) =  4.25820367510819D-02
         WT( 2, 2) =  1.15233316622473D-01
         WT( 3, 2) =  1.86800796556493D-01
         WT( 4, 2) =  2.41040339228648D-01
         WT( 5, 2) =  2.72849801912558D-01
         WT( 6, 2) =  2.82987417857491D-01
      ELSE IF (N .EQ. 6) THEN
C            (N,2N+1)=    6   13
         PT( 1) = -9.88703202612676D-01
         PT( 2) = -9.32469514203099D-01
         PT( 3) = -8.21373340865030D-01
         PT( 4) = -6.61209386466240D-01
         PT( 5) = -4.63118212475301D-01
         PT( 6) = -2.38619186083182D-01
         PT( 7) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  1.71324492379181D-01
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  3.60761573048112D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  4.67913934572650D-01
         WT( 7, 1) =                 0.D0
         WT( 1, 2) =  3.03961541198198D-02
         WT( 2, 2) =  8.36944404469064D-02
         WT( 3, 2) =  1.37320604634447D-01
         WT( 4, 2) =  1.81071994323138D-01
         WT( 5, 2) =  2.13209652271962D-01
         WT( 6, 2) =  2.33770864116995D-01
         WT( 7, 2) =  2.41072580173465D-01
      ELSE IF (N .EQ. 7) THEN
C            (N,2N+1)=    7   15
         PT( 1) = -9.91455371120814D-01
         PT( 2) = -9.49107912342697D-01
         PT( 3) = -8.64864423359769D-01
         PT( 4) = -7.41531185599360D-01
         PT( 5) = -5.86087235467687D-01
         PT( 6) = -4.05845151377365D-01
         PT( 7) = -2.07784955007895D-01
         PT( 8) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  1.29484966168870D-01
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  2.79705391489244D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  3.81830050505098D-01
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  4.17959183673426D-01
         WT( 1, 2) =  2.29353220105292D-02
         WT( 2, 2) =  6.30920926299785D-02
         WT( 3, 2) =  1.04790010322250D-01
         WT( 4, 2) =  1.40653259715526D-01
         WT( 5, 2) =  1.69004726639268D-01
         WT( 6, 2) =  1.90350578064786D-01
         WT( 7, 2) =  2.04432940075299D-01
         WT( 8, 2) =  2.09482141084727D-01
      ELSE IF (N .EQ. 9) THEN
C            (N,2N+1)=    9   19
         PT( 1) = -9.94678160677338D-01
         PT( 2) = -9.68160239507569D-01
         PT( 3) = -9.14963507249681D-01
         PT( 4) = -8.36031107326573D-01
         PT( 5) = -7.34486765183931D-01
         PT( 6) = -6.13371432700553D-01
         PT( 7) = -4.75462479112458D-01
         PT( 8) = -3.24253423403770D-01
         PT( 9) = -1.64223563614982D-01
         PT(10) =  2.96604890223891D-14
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  8.12743883615985D-02
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  1.80648160694843D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  2.60610696402914D-01
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  3.12347077039973D-01
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  3.30239355001209D-01
         WT( 1, 2) =  1.43047756438390D-02
         WT( 2, 2) =  3.96318951602612D-02
         WT( 3, 2) =  6.65181559402743D-02
         WT( 4, 2) =  9.07906816887265D-02
         WT( 5, 2) =  1.11789134684418D-01
         WT( 6, 2) =  1.30001406855341D-01
         WT( 7, 2) =  1.45239588384366D-01
         WT( 8, 2) =  1.56413527788484D-01
         WT( 9, 2) =  1.62862827440115D-01
         WT(10, 2) =  1.64896012828350D-01
      ELSE IF (N .EQ. 10) THEN
C            (N,2N+1)=   10   21
         PT( 1) = -9.95657163025811D-01
         PT( 2) = -9.73906528517077D-01
         PT( 3) = -9.30157491355708D-01
         PT( 4) = -8.65063366688911D-01
         PT( 5) = -7.80817726586413D-01
         PT( 6) = -6.79409568298979D-01
         PT( 7) = -5.62757134668601D-01
         PT( 8) = -4.33395394129214D-01
         PT( 9) = -2.94392862701457D-01
         PT(10) = -1.48874338981601D-01
         PT(11) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  6.66713443087090D-02
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  1.49451349150563D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  2.19086362515949D-01
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  2.69266719309986D-01
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  2.95524224714692D-01
         WT(11, 1) =                 0.D0
         WT( 1, 2) =  1.16946388673718D-02
         WT( 2, 2) =  3.25581623079647D-02
         WT( 3, 2) =  5.47558965743520D-02
         WT( 4, 2) =  7.50396748109199D-02
         WT( 5, 2) =  9.31254545836975D-02
         WT( 6, 2) =  1.09387158802297D-01
         WT( 7, 2) =  1.23491976262066D-01
         WT( 8, 2) =  1.34709217311474D-01
         WT( 9, 2) =  1.42775938577060D-01
         WT(10, 2) =  1.47739104901339D-01
         WT(11, 2) =  1.49445554002917D-01
      ELSE IF (N .EQ. 12) THEN
C            (N,2N+1)=   12   25
         PT( 1) = -9.96933922529593D-01
         PT( 2) = -9.81560634246620D-01
         PT( 3) = -9.50537795943120D-01
         PT( 4) = -9.04117256370380D-01
         PT( 5) = -8.43558124161156D-01
         PT( 6) = -7.69902674194224D-01
         PT( 7) = -6.84059895470057D-01
         PT( 8) = -5.87317954286593D-01
         PT( 9) = -4.81339450478153D-01
         PT(10) = -3.67831498998131D-01
         PT(11) = -2.48505748320468D-01
         PT(12) = -1.25233408511434D-01
         PT(13) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  4.71753363865584D-02
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  1.06939325995315D-01
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  1.60078328543288D-01
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  2.03167426723049D-01
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  2.33492536538344D-01
         WT(11, 1) =                 0.D0
         WT(12, 1) =  2.49147045813361D-01
         WT(13, 1) =                 0.D0
         WT( 1, 2) =  8.25771143316839D-03
         WT( 2, 2) =  2.30360840389822D-02
         WT( 3, 2) =  3.89152304692995D-02
         WT( 4, 2) =  5.36970176077562D-02
         WT( 5, 2) =  6.72509070508398D-02
         WT( 6, 2) =  7.99202753336017D-02
         WT( 7, 2) =  9.15494682950491D-02
         WT( 8, 2) =  1.01649732279060D-01
         WT( 9, 2) =  1.10022604977644D-01
         WT(10, 2) =  1.16712053501757D-01
         WT(11, 2) =  1.21626303523948D-01
         WT(12, 2) =  1.24584164536156D-01
         WT(13, 2) =  1.25556893905475D-01
      ELSE IF (N .EQ. 15) THEN
C            (N,2N+1)=   15   31
         PT( 1) = -9.98002298693400D-01
         PT( 2) = -9.87992518020345D-01
         PT( 3) = -9.67739075679141D-01
         PT( 4) = -9.37273392400584D-01
         PT( 5) = -8.97264532344082D-01
         PT( 6) = -8.48206583410331D-01
         PT( 7) = -7.90418501442467D-01
         PT( 8) = -7.24417731360095D-01
         PT( 9) = -6.50996741297419D-01
         PT(10) = -5.70972172608489D-01
         PT(11) = -4.85081863640239D-01
         PT(12) = -3.94151347077514D-01
         PT(13) = -2.99180007153169D-01
         PT(14) = -2.01194093997383D-01
         PT(15) = -1.01142066918719D-01
         PT(16) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  3.07532419961866D-02
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  7.03660474880725D-02
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  1.07159220467161D-01
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  1.39570677926137D-01
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  1.66269205816981D-01
         WT(11, 1) =                 0.D0
         WT(12, 1) =  1.86161000015526D-01
         WT(13, 1) =                 0.D0
         WT(14, 1) =  1.98431485327072D-01
         WT(15, 1) =                 0.D0
         WT(16, 1) =  2.02578241925522D-01
         WT( 1, 2) =  5.37747987292339D-03
         WT( 2, 2) =  1.50079473293162D-02
         WT( 3, 2) =  2.54608473267154D-02
         WT( 4, 2) =  3.53463607913758D-02
         WT( 5, 2) =  4.45897513247648D-02
         WT( 6, 2) =  5.34815246909279D-02
         WT( 7, 2) =  6.20095678006707D-02
         WT( 8, 2) =  6.98541213187283D-02
         WT( 9, 2) =  7.68496807577206D-02
         WT(10, 2) =  8.30805028231332D-02
         WT(11, 2) =  8.85644430562116D-02
         WT(12, 2) =  9.31265981708256D-02
         WT(13, 2) =  9.66427269836236D-02
         WT(14, 2) =  9.91735987217921D-02
         WT(15, 2) =  1.00769845523875D-01
         WT(16, 2) =  1.01330007014792D-01
      ELSE IF (N .EQ. 20) THEN
C            (N,2N+1)=   20   41
         PT( 1) = -9.98859031588275D-01
         PT( 2) = -9.93128599185077D-01
         PT( 3) = -9.81507877450248D-01
         PT( 4) = -9.63971927277896D-01
         PT( 5) = -9.40822633831758D-01
         PT( 6) = -9.12234428251306D-01
         PT( 7) = -8.78276811252285D-01
         PT( 8) = -8.39116971822200D-01
         PT( 9) = -7.95041428837557D-01
         PT(10) = -7.46331906460135D-01
         PT(11) = -6.93237656334752D-01
         PT(12) = -6.36053680726494D-01
         PT(13) = -5.75140446819709D-01
         PT(14) = -5.10867001950807D-01
         PT(15) = -4.43593175238725D-01
         PT(16) = -3.73706088715412D-01
         PT(17) = -3.01627868114913D-01
         PT(18) = -2.27785851141632D-01
         PT(19) = -1.52605465240924D-01
         PT(20) = -7.65265211334905D-02
         PT(21) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  1.76140071391723D-02
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  4.06014298003696D-02
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  6.26720483341199D-02
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  8.32767415767082D-02
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  1.01930119817242D-01
         WT(11, 1) =                 0.D0
         WT(12, 1) =  1.18194531961518D-01
         WT(13, 1) =                 0.D0
         WT(14, 1) =  1.31688638449168D-01
         WT(15, 1) =                 0.D0
         WT(16, 1) =  1.42096109318361D-01
         WT(17, 1) =                 0.D0
         WT(18, 1) =  1.49172986472609D-01
         WT(19, 1) =                 0.D0
         WT(20, 1) =  1.52753387130707D-01
         WT(21, 1) =                 0.D0
         WT( 1, 2) =  3.07358371852059D-03
         WT( 2, 2) =  8.60026985564299D-03
         WT( 3, 2) =  1.46261692569712D-02
         WT( 4, 2) =  2.03883734612667D-02
         WT( 5, 2) =  2.58821336049512D-02
         WT( 6, 2) =  3.12873067770327D-02
         WT( 7, 2) =  3.66001697582008D-02
         WT( 8, 2) =  4.16688733279735D-02
         WT( 9, 2) =  4.64348218674977D-02
         WT(10, 2) =  5.09445739237286D-02
         WT(11, 2) =  5.51951053482860D-02
         WT(12, 2) =  5.91114008806397D-02
         WT(13, 2) =  6.26532375547812D-02
         WT(14, 2) =  6.58345971336183D-02
         WT(15, 2) =  6.86486729285214D-02
         WT(16, 2) =  7.10544235534440D-02
         WT(17, 2) =  7.30306903327866D-02
         WT(18, 2) =  7.45828754004991D-02
         WT(19, 2) =  7.57044976845567D-02
         WT(20, 2) =  7.63778676720808D-02
         WT(21, 2) =  7.66007119179997D-02
      ELSE IF (N .EQ. 30) THEN
C            (N,2N+1)=   30   61
         PT( 1) = -9.99484410050492D-01
         PT( 2) = -9.96893484074594D-01
         PT( 3) = -9.91630996870406D-01
         PT( 4) = -9.83668123279688D-01
         PT( 5) = -9.73116322501127D-01
         PT( 6) = -9.60021864968258D-01
         PT( 7) = -9.44374444748561D-01
         PT( 8) = -9.26200047429226D-01
         PT( 9) = -9.05573307699910D-01
         PT(10) = -8.82560535792024D-01
         PT(11) = -8.57205233546061D-01
         PT(12) = -8.29565762382735D-01
         PT(13) = -7.99727835821841D-01
         PT(14) = -7.67777432104801D-01
         PT(15) = -7.33790062453224D-01
         PT(16) = -6.97850494793279D-01
         PT(17) = -6.60061064126623D-01
         PT(18) = -6.20526182989213D-01
         PT(19) = -5.79345235826359D-01
         PT(20) = -5.36624148141996D-01
         PT(21) = -4.92480467861775D-01
         PT(22) = -4.47033769538073D-01
         PT(23) = -4.00401254830392D-01
         PT(24) = -3.52704725530867D-01
         PT(25) = -3.04073202273621D-01
         PT(26) = -2.54636926167882D-01
         PT(27) = -2.04525116682305D-01
         PT(28) = -1.53869913608576D-01
         PT(29) = -1.02806937966733D-01
         PT(30) = -5.14718425552887D-02
         PT(31) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  7.96819249620989D-03
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  1.84664683111015D-02
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  2.87847078833273D-02
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  3.87991925696209D-02
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  4.84026728305875D-02
         WT(11, 1) =                 0.D0
         WT(12, 1) =  5.74931562175982D-02
         WT(13, 1) =                 0.D0
         WT(14, 1) =  6.59742298821788D-02
         WT(15, 1) =                 0.D0
         WT(16, 1) =  7.37559747377015D-02
         WT(17, 1) =                 0.D0
         WT(18, 1) =  8.07558952294163D-02
         WT(19, 1) =                 0.D0
         WT(20, 1) =  8.68997872010757D-02
         WT(21, 1) =                 0.D0
         WT(22, 1) =  9.21225222377777D-02
         WT(23, 1) =                 0.D0
         WT(24, 1) =  9.63687371746342D-02
         WT(25, 1) =                 0.D0
         WT(26, 1) =  9.95934205867908D-02
         WT(27, 1) =                 0.D0
         WT(28, 1) =  1.01762389748416D-01
         WT(29, 1) =                 0.D0
         WT(30, 1) =  1.02852652893556D-01
         WT(31, 1) =                 0.D0
         WT( 1, 2) =  1.38901369867700D-03
         WT( 2, 2) =  3.89046112709980D-03
         WT( 3, 2) =  6.63070391593121D-03
         WT( 4, 2) =  9.27327965951780D-03
         WT( 5, 2) =  1.18230152534964D-02
         WT( 6, 2) =  1.43697295070458D-02
         WT( 7, 2) =  1.69208891890532D-02
         WT( 8, 2) =  1.94141411939424D-02
         WT( 9, 2) =  2.18280358216092D-02
         WT(10, 2) =  2.41911620780806D-02
         WT(11, 2) =  2.65099548823332D-02
         WT(12, 2) =  2.87540487650412D-02
         WT(13, 2) =  3.09072575623878D-02
         WT(14, 2) =  3.29814470574838D-02
         WT(15, 2) =  3.49793380280601D-02
         WT(16, 2) =  3.68823646518213D-02
         WT(17, 2) =  3.86789456247276D-02
         WT(18, 2) =  4.03745389515360D-02
         WT(19, 2) =  4.19698102151642D-02
         WT(20, 2) =  4.34525397013561D-02
         WT(21, 2) =  4.48148001331625D-02
         WT(22, 2) =  4.60592382710070D-02
         WT(23, 2) =  4.71855465692992D-02
         WT(24, 2) =  4.81858617570872D-02
         WT(25, 2) =  4.90554345550298D-02
         WT(26, 2) =  4.97956834270743D-02
         WT(27, 2) =  5.04059214027823D-02
         WT(28, 2) =  5.08817958987495D-02
         WT(29, 2) =  5.12215478492588D-02
         WT(30, 2) =  5.14261285374591D-02
         WT(31, 2) =  5.14947294294517D-02
      ELSE IF (N .EQ. 40) THEN
C            (N,2N+1)=   40   81
         PT( 1) = -9.99707559258702D-01
         PT( 2) = -9.98237709710512D-01
         PT( 3) = -9.95250573446071D-01
         PT( 4) = -9.90726238699406D-01
         PT( 5) = -9.84722839864247D-01
         PT( 6) = -9.77259949983733D-01
         PT( 7) = -9.68323126854152D-01
         PT( 8) = -9.57916819213754D-01
         PT( 9) = -9.46071837162499D-01
         PT(10) = -9.32812808278623D-01
         PT(11) = -9.18149543072900D-01
         PT(12) = -9.02098806968837D-01
         PT(13) = -8.84692008701087D-01
         PT(14) = -8.65959503212213D-01
         PT(15) = -8.45923985587312D-01
         PT(16) = -8.24612230833264D-01
         PT(17) = -8.02060566140248D-01
         PT(18) = -7.78305651426486D-01
         PT(19) = -7.53379803438939D-01
         PT(20) = -7.27318255189900D-01
         PT(21) = -7.00162977487331D-01
         PT(22) = -6.71956684614155D-01
         PT(23) = -6.42739524305576D-01
         PT(24) = -6.12553889667957D-01
         PT(25) = -5.81447065829131D-01
         PT(26) = -5.49467125095116D-01
         PT(27) = -5.16660607386385D-01
         PT(28) = -4.83075801686166D-01
         PT(29) = -4.48764513638160D-01
         PT(30) = -4.13779204371586D-01
         PT(31) = -3.78171435473590D-01
         PT(32) = -3.41994090825745D-01
         PT(33) = -3.05302441735243D-01
         PT(34) = -2.68152185007231D-01
         PT(35) = -2.30598521880715D-01
         PT(36) = -1.92697580701347D-01
         PT(37) = -1.54506879379390D-01
         PT(38) = -1.16084070675229D-01
         PT(39) = -7.74865883312827D-02
         PT(40) = -3.87724175060129D-02
         PT(41) =                 0.D0
         WT( 1, 1) =                 0.D0
         WT( 2, 1) =  4.52127709857483D-03
         WT( 3, 1) =                 0.D0
         WT( 4, 1) =  1.04982845311367D-02
         WT( 5, 1) =                 0.D0
         WT( 6, 1) =  1.64210583818923D-02
         WT( 7, 1) =                 0.D0
         WT( 8, 1) =  2.22458491942125D-02
         WT( 9, 1) =                 0.D0
         WT(10, 1) =  2.79370069799963D-02
         WT(11, 1) =                 0.D0
         WT(12, 1) =  3.34601952825375D-02
         WT(13, 1) =                 0.D0
         WT(14, 1) =  3.87821679744755D-02
         WT(15, 1) =                 0.D0
         WT(16, 1) =  4.38709081856796D-02
         WT(17, 1) =                 0.D0
         WT(18, 1) =  4.86958076350388D-02
         WT(19, 1) =                 0.D0
         WT(20, 1) =  5.32278469839487D-02
         WT(21, 1) =                 0.D0
         WT(22, 1) =  5.74397690994031D-02
         WT(23, 1) =                 0.D0
         WT(24, 1) =  6.13062424929203D-02
         WT(25, 1) =                 0.D0
         WT(26, 1) =  6.48040134565959D-02
         WT(27, 1) =                 0.D0
         WT(28, 1) =  6.79120458152260D-02
         WT(29, 1) =                 0.D0
         WT(30, 1) =  7.06116473913001D-02
         WT(31, 1) =                 0.D0
         WT(32, 1) =  7.28865823957805D-02
         WT(33, 1) =                 0.D0
         WT(34, 1) =  7.47231690579957D-02
         WT(35, 1) =                 0.D0
         WT(36, 1) =  7.61103619006116D-02
         WT(37, 1) =                 0.D0
         WT(38, 1) =  7.70398181642697D-02
         WT(39, 1) =                 0.D0
         WT(40, 1) =  7.75059479784290D-02
         WT(41, 1) =                 0.D0
         WT( 1, 2) =  7.87863323894401D-04
         WT( 2, 2) =  2.20748573572679D-03
         WT( 3, 2) =  3.76522867934199D-03
         WT( 4, 2) =  5.27194271488540D-03
         WT( 5, 2) =  6.73181348520741D-03
         WT( 6, 2) =  8.19757638675139D-03
         WT( 7, 2) =  9.67540148401719D-03
         WT( 8, 2) =  1.11313216640276D-02
         WT( 9, 2) =  1.25543847685172D-02
         WT(10, 2) =  1.39625598669806D-02
         WT(11, 2) =  1.53613263591024D-02
         WT(12, 2) =  1.67345324750026D-02
         WT(13, 2) =  1.80738684088182D-02
         WT(14, 2) =  1.93876458943179D-02
         WT(15, 2) =  2.06790432735282D-02
         WT(16, 2) =  2.19381873358330D-02
         WT(17, 2) =  2.31589310133770D-02
         WT(18, 2) =  2.43456901822734D-02
         WT(19, 2) =  2.55002176031301D-02
         WT(20, 2) =  2.66157374990246D-02
         WT(21, 2) =  2.76876261110610D-02
         WT(22, 2) =  2.87183868410922D-02
         WT(23, 2) =  2.97089272777766D-02
         WT(24, 2) =  3.06543608914116D-02
         WT(25, 2) =  3.15512236191153D-02
         WT(26, 2) =  3.24009825076059D-02
         WT(27, 2) =  3.32040443412576D-02
         WT(28, 2) =  3.39568628342097D-02
         WT(29, 2) =  3.46569358434976D-02
         WT(30, 2) =  3.53051447086219D-02
         WT(31, 2) =  3.59016027836283D-02
         WT(32, 2) =  3.64438265303411D-02
         WT(33, 2) =  3.69301695340487D-02
         WT(34, 2) =  3.73611800254692D-02
         WT(35, 2) =  3.77368012630936D-02
         WT(36, 2) =  3.80554637788524D-02
         WT(37, 2) =  3.83163240051747D-02
         WT(38, 2) =  3.85197417499508D-02
         WT(39, 2) =  3.86655554391411D-02
         WT(40, 2) =  3.87530293787524D-02
         WT(41, 2) =  3.87821047642829D-02
      ELSE
         WRITE (6, 600) N
         STOP 'KG1 1'
      END IF
      MP = M + 1
      DO 20 I = 1,N
         II = MP - I
         PT(II) = -PT(I)
         DO 10 J = 1,2
            WT(II,J) = WT(I,J)
   10    CONTINUE
   20 CONTINUE
      RETURN
  600 FORMAT (2X,T5,'Error: Number of quadrature points is', I3/ 12X,
     * ' choose a value from 2-7,9,10,12,15,20,30,40 and try again')
      END
