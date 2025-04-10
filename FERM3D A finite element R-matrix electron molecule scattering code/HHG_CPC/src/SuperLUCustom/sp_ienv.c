/**************************************************************
 * 文件/File: sp_ienv.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:07
 **************************************************************
 */

/**************************************************************
 * 文件/File: sp_ienv.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:43
 **************************************************************
 */

/*
 * File name:		sp_ienv.c
 * History:             Modified from lapack routine ILAENV
 */
int
sp_ienv(int ispec)
{
/*
    Purpose   
    =======   

    sp_ienv() is inquired to choose machine-dependent parameters for the
    local environment. See ISPEC for a description of the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    Arguments   
    =========   

    ISPEC   (input) int
            Specifies the parameter to be returned as the value of SP_IENV.   
            = 1: the panel size w; a panel consists of w consecutive
	         columns of matrix A in the process of Gaussian elimination.
		 The best value depends on machine's cache characters.
            = 2: the relaxation parameter relax; if the number of
	         nodes (columns) in a subtree of the elimination tree is less
		 than relax, this subtree is considered as one supernode,
		 regardless of the their row structures.
            = 3: the maximum size for a supernode;
	    = 4: the minimum row dimension for 2-D blocking to be used;
	    = 5: the minimum column dimension for 2-D blocking to be used;
	    = 6: the estimated fills factor for L and U, compared with A;
	    
   (SP_IENV) (output) int
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if SP_IENV = -k, the k-th argument had an illegal value. 
  
    ===================================================================== 
*/
    int i;

    switch (ispec) {
	case 1: return (8);
	case 2: return (5);
	case 3: return (100);
	case 4: return (400);
	case 5: return (40);
        case 6: return (16);
    }

    /* Invalid value for ISPEC */
    i = 1;
    xerbla_("sp_ienv", &i);
    return 0;

} /* sp_ienv_ */

