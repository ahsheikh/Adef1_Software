//  Here we define user-defined functions
#include "MGSLP.h"
#include "petscksp.h"
#include "petscpcmg.h"
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGSLPShellPCCreate"

PetscErrorCode MGSLPShellPCCreate(MGSLPShellPC **shell)
{
   MGSLPShellPC *newctx;
   PetscErrorCode ierr;
      ierr = PetscNew(MGSLPShellPC,&newctx); CHKERRQ(ierr);
   *shell = newctx;
   return 0;
}
 /* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGSLPShellPCSetUp"
 PetscErrorCode MGSLPShellPCSetUp(PC pc,GridCtx* grid,Vec x)
 {
    PetscErrorCode ierr;
    Vec xforM, x_pcshift;
    int k,petsc_level,  levels,levelsforM,sizeA;
    KSP coarsegridksp,ksp1;
    PC  coarsegridpc,pc1; 
    pcshift	  *shellpcshift;
    MGforM	  *shellforM;
    Vec xPC;

    ierr = PCSetType(pc,PCMG);CHKERRQ(ierr);
    levels = grid[0].lev+1;
    ierr = PetscOptionsGetInt(PETSC_NULL, "-levels", &levels, 0);CHKERRQ(ierr);
    ierr = PCMGSetLevels(pc,levels,PETSC_NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
    ierr = PCMGSetCycleType(pc, PC_MG_CYCLE_V); CHKERRQ(ierr);
//     ierr = PCMGMultiplicativeSetCycles(pc, 1);CHKERRQ(ierr);

for (k=1; k<levels; k++){  
     petsc_level = levels - k;  /*petsc_level is current level.dont give coarsest = 0.*/
     levelsforM = petsc_level+1; 
//     ....Get pre-smoothing context....
    ierr = PCMGGetSmootherDown(pc,petsc_level,&grid[k].ksp_pre); CHKERRQ(ierr); 
    ierr = PCMGGetSmootherUp(pc,petsc_level,&grid[k].ksp_post); CHKERRQ(ierr); 
    
    ierr = MatGetSize(grid[k+1].A, &sizeA, &sizeA); CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_WORLD,&grid[k].x); CHKERRQ(ierr);
    ierr = VecSetType(grid[k].x, VECSEQ); CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_WORLD,&grid[k].b); CHKERRQ(ierr);
    ierr = VecSetType(grid[k].b, VECSEQ); CHKERRQ(ierr); 
    ierr = VecCreate(MPI_COMM_WORLD,&grid[k].r); CHKERRQ(ierr);
    ierr = VecSetType(grid[k].r,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(grid[k].x,PETSC_DECIDE,sizeA); CHKERRQ(ierr);
    ierr = VecSetSizes(grid[k].r,PETSC_DECIDE,sizeA);  CHKERRQ(ierr);
    ierr = VecSetSizes(grid[k].b,PETSC_DECIDE,sizeA);  CHKERRQ(ierr);
//     ...set ksp_pre context....
    ierr = KSPSetOperators(grid[k].ksp_pre,grid[k].A,grid[k].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(grid[k].ksp_pre,&grid[k].pc_pre);CHKERRQ(ierr);
    ierr = KSPSetType(grid[k].ksp_pre, KSPRICHARDSON);CHKERRQ(ierr);
//     ierr = KSPSetTolerances(grid[k].ksp_pre, 1e-12, 1e-50, 1e7,0);CHKERRQ(ierr); 
    ierr = PCSetType(grid[k].pc_pre, PCNONE);CHKERRQ(ierr);
      //     ========================================================================|
//      ierr = PCSetType(grid[k].pc_pre,PCSHELL); CHKERRQ(ierr); 			   //|
//      ierr = MGforMCreate(&shellforM); CHKERRQ(ierr);				   //|
//      ierr = PCShellSetApply(grid[k].pc_pre,MGforMApply); CHKERRQ(ierr);		   //|
//      ierr = PCShellSetContext(grid[k].pc_pre,shellforM); CHKERRQ(ierr);		   //|
//      ierr = PCShellSetDestroy(grid[k].pc_pre,MGforMDestroy);CHKERRQ(ierr);	   //|
//      ierr = MGforMSetUp(grid[k].pc_pre,&grid[0],xforM,levelsforM ); CHKERRQ(ierr); //|
     //     =========================================================================|
    
    ierr = PCMGSetNumberSmoothDown(pc,0); CHKERRQ(ierr);
// ierr = KSPView(grid[k].ksp_prevery,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
//     ....set KSP_POST context....
    ierr = KSPSetOperators(grid[k].ksp_post,grid[k].A,grid[k].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(grid[k].ksp_post,&grid[k].pc_post); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(grid[k].ksp_post, PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetType(grid[k].ksp_post,KSPRICHARDSON); CHKERRQ(ierr);
//     ierr = KSPSetTolerances(grid[k].ksp_post,1e-12, 1e-50, 1e7,1); CHKERRQ(ierr);  
    ierr = PCSetType(grid[k].pc_post, PCLU); CHKERRQ(ierr); 
    ierr = PCMGSetNumberSmoothUp(pc,1); CHKERRQ(ierr);
    //     =========================================================
//      ierr = PCSetType(grid[k].pc_post,PCSHELL); CHKERRQ(ierr); 
//      ierr = pcshiftCreate(&shellpcshift); CHKERRQ(ierr);
//      ierr = PCShellSetApply(grid[k].pc_post,pcshiftApply); CHKERRQ(ierr);
//      ierr = PCShellSetContext(grid[k].pc_post,shellpcshift); CHKERRQ(ierr);
//      ierr = PCShellSetDestroy(grid[k].pc_post,pcshiftDestroy);CHKERRQ(ierr);
//      ierr = PCShellSetName(grid[k].pc_post,"MyPreconditioner");CHKERRQ(ierr);  
//      ierr = pcshiftSetUp(grid[k].pc_post,&grid[0],xPC); CHKERRQ(ierr);
//      =========================================================
// ierr = KSPView(grid[k].ksp_post, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
    //     ==============================
    ierr = PCMGSetX(pc,petsc_level-1,grid[k].x);CHKERRQ(ierr); 
    ierr = PCMGSetRhs(pc,petsc_level-1,grid[k].b);CHKERRQ(ierr); 
//     ierr = PCMGSetR(pc,petsc_level,grid[k].r);CHKERRQ(ierr);
    ierr = PCMGSetResidual(pc,petsc_level,PCMGDefaultResidual,grid[k].A); CHKERRQ(ierr);
//     ....Create interpolation between the levels....
    ierr = PCMGSetInterpolation(pc,petsc_level,grid[k].P);CHKERRQ(ierr);
    ierr = PCMGSetRestriction(pc,petsc_level,grid[k].R);CHKERRQ(ierr);  
}
  /*....Set coarse grid solver....*/  
    ierr = PCMGGetCoarseSolve(pc,&coarsegridksp); CHKERRQ(ierr); 
    ierr = KSPSetFromOptions(coarsegridksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(coarsegridksp,grid[levels].A,grid[levels].A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(coarsegridksp,&coarsegridpc); CHKERRQ(ierr);
    ierr = KSPSetType(coarsegridksp, KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPSetTolerances(coarsegridksp, 1e-12, 1e-50, 1e7,1);CHKERRQ(ierr); 
    ierr = PCSetType(coarsegridpc, PCLU); CHKERRQ(ierr); 
   //   ierr = KSPView(coarsegridksp, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
   //   ierr = KSPSetMonitor(coarsegridksp,KSPDefaultMonitor,PETSC_NULL, 
   //          PETSC_NULL); CHKERRQ(ierr); 
 /*..Allow above criterea to be overwritten..*/
    ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
//     ====================================================
//     OR SET BELOW THE PETSC ILU PRECONDITIONER 
//     ====================================================
// ierr = PCSetType(pc, PCILU); CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_SELF, "Inside user defined preconditioner\n");
ierr = PCSetUp(pc); CHKERRQ(ierr);
 return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGSLPShellPCApply"

PetscErrorCode MGSLPShellPCApply(PC pc,Vec x,Vec y)
 {
      PetscErrorCode ierr;

ierr = PCApply(pc, x , y); CHKERRQ(ierr);
     
   return 0;
 }
 
#undef __FUNC__
#define __FUNC__ "MGSLPShellPCDestroy"
 PetscErrorCode MGSLPShellPCDestroy(PC pc) 
{

  /*..Free PCShell context..*/
//   PetscFree(shell);

  return 0;
}
// ===================================================================================================================
//                         here we start to defined MG approximation for M 
// ===================================================================================================================

/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGSforMCreate"

/*.. MGSLPShellPCCreate - This routine creates a user-defined
     preconditioner context.

     Output Parameter:
     shell - user-defined preconditioner context..*/

PetscErrorCode MGforMCreate(MGforM **shellforM)
{
   MGforM *newctxforM; 
   PetscErrorCode ierr;
 
     ierr = PetscNew(MGforM,&newctxforM); CHKERRQ(ierr);
   *shellforM = newctxforM; 
   return 0; 
}
  
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGforMSetUp"
 PetscErrorCode MGforMSetUp(PC pcforM,GridCtx* grid,Vec x,int levelsforM)
 {
    PetscErrorCode ierr;
// //     Mat MHtemp;
//     int myind, ind,sizeM, petsc_levelforM, levels; //petsc_level, petsc_levelForM, levels,levelsForM;
//     PetscReal    omega;
//     KSP  CGkspForM; 
//     PC   CGpcForM; 
//      omega = .7;               /*Relaxation factor*/
// // PetscPrintf(PETSC_COMM_WORLD,"Number of INNER LEVELS for M = %d\n",levelsforM);
// 
    ierr = PCSetType(pcforM,PCLU);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_SELF, "Here MGforM applies PCLU on M \n");
// // ================================PCMG  structure for M ===============================================
//     
//     ierr = PCSetType(pcforM,PCMG);CHKERRQ(ierr);
//     ierr = PetscOptionsGetInt(PETSC_NULL, "-levelsforM", &levelsforM, 0);CHKERRQ(ierr);
//     ierr = PCMGSetLevels(pcforM,levelsforM,PETSC_NULL);CHKERRQ(ierr);
//     ierr = PCMGSetType(pcforM, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
//     ierr = PCMGMultiplicativeSetCycles(pcforM,1);CHKERRQ(ierr);
//     ierr = PCMGSetNumberSmoothDown(pcforM,1); CHKERRQ(ierr);
//     ierr = PCMGSetNumberSmoothUp(pcforM,1); CHKERRQ(ierr);
//     levels = grid[0].lev+1;
// for (myind=1; myind<levelsforM; myind++){  
//      petsc_levelforM = levelsforM - myind;  /*petsc_level is current level.dont give coarsest = 0.*/
// ind = levels - petsc_levelforM; 
// //     ....Get pre-smoothing context....
//     ierr = PCMGGetSmootherDown(pcforM,petsc_levelforM,&grid[ind].ksp_preForM); CHKERRQ(ierr); 
// //     ierr = KSPView(grid[ind].ksp_preForM,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
//     ierr = PCMGGetSmootherUp(pcforM,petsc_levelforM,&grid[ind].ksp_postForM); CHKERRQ(ierr); 
// //    
//     ierr = MatGetSize(grid[ind+1].M, &sizeM, &sizeM); CHKERRQ(ierr);
//     ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].x_M); CHKERRQ(ierr);
//     ierr = VecSetType(grid[ind].x_M, VECSEQ); CHKERRQ(ierr);
//     ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].b_M); CHKERRQ(ierr);
//     ierr = VecSetType(grid[ind].b_M, VECSEQ); CHKERRQ(ierr); 
//     ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].r_M); CHKERRQ(ierr);
//     ierr = VecSetType(grid[ind].r_M,VECSEQ); CHKERRQ(ierr);
//     ierr = VecSetSizes(grid[ind].x_M,PETSC_DECIDE,sizeM); CHKERRQ(ierr);
//     ierr = VecSetSizes(grid[ind].r_M,PETSC_DECIDE,sizeM);  CHKERRQ(ierr);
//     ierr = VecSetSizes(grid[ind].b_M,PETSC_DECIDE,sizeM);  CHKERRQ(ierr);
// // //     ierr = VecDuplicate(grid[k].x,&grid[k].b); CHKERRQ(ierr); 
// // //     ierr = VecDuplicate(grid[k].x,&grid[k].r); CHKERRQ(ierr);
// //   
// //     ...set ksp_preForM context....
//     ierr = KSPSetOperators(grid[ind].ksp_preForM,grid[ind].M,grid[ind].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//     ierr = KSPGetPC(grid[ind].ksp_preForM,&grid[ind].pc_preForM);CHKERRQ(ierr);
//     ierr = KSPSetType(grid[ind].ksp_preForM, KSPRICHARDSON);CHKERRQ(ierr);
//     ierr = KSPSetTolerances(grid[ind].ksp_preForM, 1e-12, 1e-50, 1e7,1);CHKERRQ(ierr); 
// //     PetscPrintf(PETSC_COMM_SELF, "Performing PREsmoothing with Damped Jacobi on M @ Level = %d\n",levelsforM); 
//     ierr = KSPRichardsonSetScale(grid[ind].ksp_preForM, omega); CHKERRQ(ierr);
//     ierr = PCSetType(grid[ind].pc_preForM, PCJACOBI);CHKERRQ(ierr); 
// //     ....set KSP_POST context....
//     ierr = KSPSetOperators(grid[ind].ksp_postForM,grid[ind].M,grid[ind].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//     ierr = KSPGetPC(grid[ind].ksp_postForM,&grid[ind].pc_postForM); CHKERRQ(ierr);
//     ierr = KSPSetInitialGuessNonzero(grid[ind].ksp_postForM, PETSC_TRUE);CHKERRQ(ierr);
//     ierr = KSPSetType(grid[ind].ksp_postForM, KSPRICHARDSON); CHKERRQ(ierr);
//     ierr = KSPSetTolerances(grid[ind].ksp_postForM,1e-12, 1e-50, 1e7,1); CHKERRQ(ierr);  
// //     PetscPrintf(PETSC_COMM_SELF, "Performing POSTsmoothing with Damped Jacobi on M @ Level = %d\n",levelsforM); 
//     ierr = KSPRichardsonSetScale(grid[ind].ksp_preForM, omega); CHKERRQ(ierr);
//     ierr = PCSetType(grid[ind].pc_postForM, PCJACOBI); CHKERRQ(ierr);   
// // //     ============================================
// //     
//     ierr = PCMGSetX(pcforM,petsc_levelforM-1,grid[ind].x_M);CHKERRQ(ierr); 
//     ierr = PCMGSetRhs(pcforM,petsc_levelforM-1,grid[ind].b_M);CHKERRQ(ierr); 
// //     ierr = PCMGSetR(pc,petsc_level,grid[k].r);CHKERRQ(ierr);
//     ierr = PCMGSetResidual(pcforM,petsc_levelforM,PCMGDefaultResidual,grid[ind].M); CHKERRQ(ierr);
// //     ....Create interpolation between the levels....
//     ierr = PCMGSetInterpolation(pcforM,petsc_levelforM,grid[ind].P);CHKERRQ(ierr);
//     ierr = PCMGSetRestriction(pcforM,petsc_levelforM,grid[ind].R);CHKERRQ(ierr);  
// //      ierr = PCMGSetInterpolation(pc,1,R);CHKERRQ(ierr);
// //      ierr = PCMGSetRestriction(pc,1,P);CHKERRQ(ierr);  
// }
//   /*....Set coarse grid solver....*/  
//     ierr = PCMGGetCoarseSolve(pcforM,&CGkspForM); CHKERRQ(ierr); 
//     ierr = KSPSetFromOptions(CGkspForM); CHKERRQ(ierr);
//     ierr = KSPSetOperators(CGkspForM,grid[levels].A,grid[levelsforM].A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//     ierr = KSPGetPC(CGkspForM,&CGpcForM); CHKERRQ(ierr);
//     ierr = KSPSetType(CGkspForM, KSPPREONLY); CHKERRQ(ierr);
//     ierr = PCSetType(CGpcForM, PCLU); CHKERRQ(ierr); 
//    //   ierr = KSPView(coarsegridksp, VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
//    //   ierr = KSPSetMonitor(coarsegridksp,KSPDefaultMonitor,PETSC_NULL, 
//    //          PETSC_NULL); CHKERRQ(ierr); 
//  /*..Allow above criterea to be overwritten..*/
//     ierr = PCSetFromOptions(pcforM); CHKERRQ(ierr);
// // //     ====================================================
// // //     OR SET BELOW THE PETSC ILU PRECONDITIONER 
// // //     ====================================================
// // // ierr = PCSetType(pc, PCILU); CHKERRQ(ierr);
// // PetscPrintf(PETSC_COMM_SELF, "Inside user defined preconditioner\n");
// // =====================================================================================================
    ierr = PCSetUp(pcforM); CHKERRQ(ierr);
 return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGforMApply"

PetscErrorCode MGforMApply(PC pcforM,Vec x,Vec y)
 {
      PetscErrorCode ierr;

ierr = PCApply(pcforM, x , y); CHKERRQ(ierr);
     
   return 0;
 }
 
#undef __FUNC__
#define __FUNC__ "MGforMDestroy"
 PetscErrorCode MGforMDestroy(PC pcforM) 
{

  /*..Free PCShell context..*/
//   PetscFree(shell);

  return 0;
}

// ===================================================================================================================
//                         here we start to defined MG approximation for M 
// ===================================================================================================================

/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "pcshiftCreate"

PetscErrorCode pcshiftCreate(pcshift **shellpcshift)
{
   pcshift *newctxpcshift; 
   PetscErrorCode ierr;
      ierr = PetscNew(pcshift,&newctxpcshift); CHKERRQ(ierr);
   *shellpcshift = newctxpcshift; 
   return 0; 
}
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "pcshiftSetUp"
 PetscErrorCode pcshiftSetUp(PC pcshift,GridCtx* grid,Vec x)
 {
    PetscErrorCode ierr;
   ierr = PCSetUp(pcshift); CHKERRQ(ierr);
 return 0;
}
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "pcshiftApply"
 PetscErrorCode pcshiftApply(PC pcshift,Vec x,Vec y)
 {
      PetscErrorCode ierr;
ierr = PCApply(pcshift, x , y); CHKERRQ(ierr);    
   return 0;
 }
/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGforMDestroy"
 PetscErrorCode pcshiftDestroy(PC pcshift) 
{
  /*..Free PCShell context..*/
//   PetscFree(shell);

  return 0;
}