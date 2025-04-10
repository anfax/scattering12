/**************************************************************
 * 文件/File: dSLUsolve.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:07
 **************************************************************
 */

/**************************************************************
 * 文件/File: dSLUsolve.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:43
 **************************************************************
 */

#include <stdlib.h>
#include <stdio.h>

#include "dsp_defs.h"
#include "util.h"
#include "Cnames.h"


int
dslusolve_(int *n, int *nnz, int *nrhs, double *values,
                int *rowind, int *colptr, double *b, int *ldb,
                int *info, char *save_fact, int *verbosity_in)

/* save_fact:  'Y'--> save factorization
**             'N'--> discard factorization
*/
{
    static SuperMatrix L, U;
    SuperMatrix A, B;
    static SCformat *Lstore;
    static NCformat *Ustore;
    static int      *perm_c; /* column permutations from partial pivoting */
    static int      *perm_r; /* row permutations from partial pivoting */
    static int      panel_size, permc_spec, i;
    static mem_usage_t   mem_usage;
    char trans[1];
    double *utime, t1;
    flops_t  *ops;
    extern SuperLUStat_t SuperLUStat;
    static int          Factored = 0;
    int Fortran, Factor, Solve, Free, Verbosity;
    int relax;

    *trans = 'N';
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(panel_size, relax);
    utime = SuperLUStat.utime;

    if ( *nrhs > 0 ) {
      Solve = 1;
    } else {
      Solve = 0;
    }
    if ( Factored ) {
      Factor = 0;
    } else {
      Factor = 1;
    }
    if ( Factor || Solve ) {
      if ( (*verbosity_in < 0) || (*verbosity_in > 10) ) {
        Verbosity = 1;
      } else {
        Verbosity = *verbosity_in;
      }
    } else {
      Verbosity = 0;
    }
    
    switch (*save_fact) {
    case 'Y': case 'y':
      Free = 0;
      break;
    case 'N': case 'n': case'1':
      Free = 1;
      break;
    default:
      printf("SuperLU: unknown option for save_fact\n");
      Factor = 0; Solve = 0; Free = 1;
    }

    Fortran = (colptr[0]==1);
    if (Factor&&Fortran) {
      /* Adjust to 0-based indexing */
      for (i = 0; i < *nnz; ++i) --rowind[i];
      for (i = 0; i <= *n; ++i) --colptr[i];
    }
    
    if (Factor) {
      dCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr, NC, _D, GE);
    }
    if (Solve) {
      dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, DN, _D, GE);
    } else {
      dCreate_Dense_Matrix(&B, 1, 1, b, 1, DN, _D, GE);
    }

    if (Factor) {
      if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
      if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    }
    if (Factor) {
      /*
       * Get column permutation vector perm_c[], according to permc_spec:
       *   permc_spec = 0: use the natural ordering 
       *   permc_spec = 1: use minimum degree ordering on structure of A'*A
       *   permc_spec = 2: use minimum degree ordering on structure of A'+A
       */
      t1 = SuperLU_timer_();
      permc_spec = 2;
      get_perm_c(permc_spec, &A, perm_c);
      utime[COLPERM] = SuperLU_timer_() - t1;
    }
    *info = 0;
    
    if (Factor) {
      /*dgssv(&A, perm_c, perm_r, &L, &U, &B, info, Solve);*/
/*** this "if" block needs to be checked with zgssv/dgssv ***/
      char refact[1];
      int lwork = 0, *etree;
      SuperMatrix AC;
      double diag_pivot_thresh = 0.0;
      double drop_tol = 0;
      DNformat *Bstore;

      *refact = 'N';
      Bstore = B.Store;

      if ( A.nrow != A.ncol || A.nrow < 0 ||
           A.Stype != NC || A.Dtype != _D || A.Mtype != GE )
          *info = -1;

      if ( *info != 0 ) {
          i = -(*info);
          xerbla_("dgssv", &i);
          return i;
      }

      if ( !(etree = intMalloc(A.ncol)) ) ABORT("Malloc fails for etree[].");

      t1 = SuperLU_timer_();
      sp_preorder(refact, &A, perm_c, etree, &AC);
      utime[ETREE] = SuperLU_timer_() - t1;

      t1 = SuperLU_timer_();
      dgstrf(refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
             etree, NULL, lwork, perm_r, perm_c, &L, &U, info);
      utime[FACT] = SuperLU_timer_() - t1;

      SUPERLU_FREE (etree);
      Destroy_CompCol_Permuted(&AC);

      if ( *info == 0 ) {
        Factored = 1;
        panel_size = sp_ienv(1);

        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
        if ( Verbosity > 2 ) {
          printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
          printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
          printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);

          dQuerySpace(&L, &U, panel_size, &mem_usage);
          printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
                  mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
                  mem_usage.expansions);
        }
      } else {
          printf("dgssv() error returns INFO= %d\n", *info);
	  if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, panel_size, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
		   mem_usage.expansions);
	  }
      }
        if ( Verbosity > 1 ) {
            utime = SuperLUStat.utime;
            ops   = SuperLUStat.ops;
            printf("GetPermC time = %8.2f\n", utime[COLPERM]);
            printf("ETREE time = %8.2f\n", utime[ETREE]);
            printf("Factor time = %8.2f\n", utime[FACT]);
            if ( utime[FACT] != 0.0 )
                printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
                       ops[FACT]*1e-6/utime[FACT]);
        }
      Destroy_SuperMatrix_Store(&A);
    } 
    if ((*info == 0) && Solve) {
        if ( B.ncol < 0 ) {
            *info = -6;
            printf("problem with B.ncol\n");
        }
        if ( B.Stype != DN ) {
            *info = -6;
            printf("problem with B.Stype\n");
        }
        if ( B.Dtype != _D ) {
            *info = -6;
            printf("problem with B.Dtype\n");
        }
        if ( B.Mtype != GE ) {
            *info = -6;
            printf("problem with B.Mtype\n");
        }
        t1 = SuperLU_timer_();
	/* Solve the system A*X=B, overwriting B with X. */
	dgstrs (trans, &L, &U, perm_r, perm_c, &B, info);
        utime[SOLVE] = SuperLU_timer_() - t1;

        if ( Verbosity > 3 ) {
            utime = SuperLUStat.utime;
            ops   = SuperLUStat.ops;
            printf("Solve time = %8.2f\n", utime[SOLVE]);
            if ( utime[SOLVE] != 0.0 )
                printf("Solve flops = %e\tMflops = %8.2f\n", ops[SOLVE],
                       ops[SOLVE]*1e-6/utime[SOLVE]);
        }
    }

    if ( Free ) {
      Factored = 0;
      SUPERLU_FREE (perm_c);
      SUPERLU_FREE (perm_r);
      Destroy_SuperNode_Matrix(&L);
      Destroy_CompCol_Matrix(&U);
    }
    StatFree();
    Destroy_SuperMatrix_Store(&B);
    if (Factor&&Fortran) {
      /* Restore to 1-based indexing */
      for (i = 0; i < *nnz; ++i) ++rowind[i];
      for (i = 0; i <= *n; ++i) ++colptr[i];
    }
    return 0;
}
