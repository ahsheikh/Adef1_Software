#ifndef PETSCFUNC_H
#define PETSCFUNC_H

/* #include "sles.h" */ 
#include "petscksp.h"



// /*..Multigrid structure for PETSc..*/ 
// 
// /*....Maximum number of levels to be used in SAMG....*/ 
#define MAX_LEVELS 25 
// 
typedef struct{
  /*..Implementation notes
    1 - The menber A is not stored on level 1 (the finest level in SAMG 
        ordering) to avoid unnecessary memory useage. 
  */ 
   KSP ksp_pre, ksp_preForM;  
   KSP ksp_post, ksp_postForM;
   PC pc_pre, pc_post, pc_preForM, pc_postForM;
  Mat  A, M, P, R;
  Mat  Interp; 
  Vec  x,r,b, upd_b,y, x_M,b_M,r_M,xforM; 
  int  size, nH,lev,iter; /*..Number of variables on level..*/ 
  //  int  debug;
} GridCtx;
// ============================================================================
//        Define function which is used recursively on coarse grid matrix.  
// ============================================================================
typedef struct{
 
  Mat  A, B, C;
//   Mat  Interp; 
  Vec  x, b, upd_b, r, y, b_y, r_y; 
  int  size; /*..Number of variables on level..*/

} MGSLPShellPC; 

/*..interface to MGSLPShellPC..*/ 
extern PetscErrorCode MGSLPShellPCCreate(MGSLPShellPC **shell); 
extern PetscErrorCode MGSLPShellPCSetUp(PC pc,GridCtx* grid,Vec x, PetscInt k);	
extern PetscErrorCode MGSLPShellPCApply(PC pc, Vec r, Vec z); 
extern PetscErrorCode MGSLPShellPCDestroy(PC pc); 
// extern int RamgGetParam(RAMG_PARAM *ramg_param);

// =========================================================================================
//                        Here we head function MGforM to apporoximate M 
//  =========================================================================================

typedef struct{
 
  //   KSP ksp_pre;  
//   KSP ksp_post;
  Mat  A, B, C;
//   Mat  Interp; 
  Vec  x, b, upd_b, r, y, b_y, r_y; 
  int  size, levelsforM; /*..Number of variables on level..*/
/*    double            *A; */
    int               *IA; 
    int               *JA;
    double            *U_APPROX; 
    double            *RHS;
    int               *IG;    
    struct RAMG_PARAM *PARAM;    
} MGforM;



extern PetscErrorCode MGforMCreate(MGforM **shellforM); 
extern PetscErrorCode MGforMSetUp(PC pcforM,GridCtx* grid,Vec x,PetscInt k);	
extern PetscErrorCode MGforMApply(PC pcforM, Vec r, Vec z); 
extern PetscErrorCode MGforMDestroy(PC pcforM); 


// =========================================================================================
//                  Here we head function "pcshift" to appply * of A + shift to 1 
// =========================================================================================

// typedef struct{
 
//     KSP ksp_pre;  
//   KSP ksp_post;
//   Mat  A, B, C;
// //   Mat  Interp; 
//   Vec  x, b, upd_b, r, y, b_y, r_y; 
//   int  size, levelsforM; /*..Number of variables on level..*/
// /*    double            *A; */
//     int               *IA; 
//     int               *JA;
//     double            *U_APPROX; 
//     double            *RHS;
//     int               *IG;    
//     struct RAMG_PARAM *PARAM;    
// } pcshift;



// extern PetscErrorCode pcshiftCreate(pcshift **shellpcshift); 
// extern PetscErrorCode pcshiftSetUp(PC pcshift,GridCtx* grid,Vec x);	
// extern PetscErrorCode pcshiftApply(PC pcshift, Vec r, Vec z); 
// extern PetscErrorCode pcshiftDestroy(PC pcshift); 







// 
// /*..Structure used in the interface to SAMG..*/ 
// typedef struct{
//     double *A; 
//     int    *IA; 
//     int    *JA; 
//     struct SAMG_PARAM *PARAM; 
//     int    LEVELS;           /* Number of levels created */   
// } SamgShellPC; 
// 
// /*..Interface to SAMG..*/ 
// extern int SamgShellPCCreate(SamgShellPC **shell); 
// extern int SamgShellPCSetUp(SamgShellPC *shell, Mat pmat);
// extern int SamgShellPCApply(void *ctx, Vec r, Vec z); 
// extern int SamgShellPCDestroy(SamgShellPC *shell); 
// extern int SamgGetParam(SAMG_PARAM *samg_param);
// 
 
// 
// /*..Level 2 routine to get coarser level matrices..*/ 
//  extern int SamgGetCoarseMat(int level, int ia_shift, int ja_shift, 
//                             Mat* coarsemat, void* ctx);
// /*..Level 2 routine to get interpolation operators..*/ 
//  extern int SamgGetInterpolation(int level, int iw_shift, int jw_shift,
//                                 Mat* interpolation, void* ctx) ;
// 
// /*..Parse SAMG hierarchy to PETSc..*/ 
// extern int SamgGetGrid(int levels, int numnodes, int numnonzero, 
//                        GridCtx* grid, void* ctx); 
// /*..Check parsing..*/ 
// extern int SamgCheckGalerkin(int levels, Mat A, GridCtx* grid, 
//                       void* ctx);

#endif//PETSCFUNC_H
