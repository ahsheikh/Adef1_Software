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
 PetscErrorCode MGSLPShellPCSetUp(PC pc,GridCtx* grid,Vec x, PetscInt k)
 {
   PetscErrorCode ierr; 
   PetscBool	LUonM=PETSC_FALSE;
    Vec  xforC; /*xforM;*/
    int coarselevel,petsc_level,sizeA,levels = 2;
    KSP coarsegridksp;
    PC  coarsegridpc; 
    MGSLPShellPC *shell;
//     pcshift	  *shellpcshift;
    MGforM	  *shellforM;
//     Vec xPC;

    ierr = PCSetType(pc,PCMG);CHKERRQ(ierr);
    coarselevel = grid[0].lev;
//     coarselevel =  5;   levels = 2; 
    ierr = PetscOptionsGetInt(PETSC_NULL, "-coarselevel", &coarselevel, 0);CHKERRQ(ierr);
//     ierr = PetscOptionsGetInt(PETSC_NULL, "-levels", &levels, 0);CHKERRQ(ierr);
    ierr = PCMGSetLevels(pc,levels,PETSC_NULL);CHKERRQ(ierr);
    ierr = PCMGSetType(pc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
    ierr = PCMGSetCycleType(pc, PC_MG_CYCLE_V); CHKERRQ(ierr);
//     ierr = PCMGMultiplicativeSetCycles(pc, 1);CHKERRQ(ierr);
//     PetscPrintf(PETSC_COMM_WORLD,"\n  ------Number of LEVELS  =  %d\n",levels);
//     PetscPrintf(PETSC_COMM_WORLD,"\n  ------COARSELEVEL @  =  %d\n",coarselevel);
// for (k=1; k<levels; k++){  
//      petsc_level = levels - k;  /*petsc_level is current level.dont give coarsest = 0.*/
     petsc_level = 1; 
//      levelsforM = petsc_level+1; 
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
     //     =========================================================
    ierr = PetscOptionsGetBool(PETSC_NULL,"-LUonM",&LUonM,PETSC_NULL);CHKERRQ(ierr);
  if (LUonM) {
//     ierr = KSPSetTolerances(grid[k].ksp_post,1e-12, 1e-50, 1e7,1); CHKERRQ(ierr);  
    ierr = PCSetType(grid[k].pc_post, PCLU); CHKERRQ(ierr);
      } else {
     ierr = PCSetType(grid[k].pc_post,PCSHELL); CHKERRQ(ierr); 
     ierr = MGforMCreate(&shellforM); CHKERRQ(ierr);
     ierr = PCShellSetApply(grid[k].pc_post,MGforMApply); CHKERRQ(ierr);
     ierr = PCShellSetContext(grid[k].pc_post,shellforM); CHKERRQ(ierr);
     ierr = PCShellSetDestroy(grid[k].pc_post,MGforMDestroy);CHKERRQ(ierr);
     ierr = MGforMSetUp(grid[k].pc_post,&grid[0],grid[k].xforM,k); CHKERRQ(ierr); 
     
//      ierr = KSPView(grid[k].ksp_post, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
      }
      ierr = PCMGSetNumberSmoothUp(pc,1); CHKERRQ(ierr);
//      =========================================================
    ierr = PCMGSetX(pc,petsc_level-1,grid[k].x);CHKERRQ(ierr); 
    ierr = PCMGSetRhs(pc,petsc_level-1,grid[k].b);CHKERRQ(ierr); 
//     ierr = PCMGSetR(pc,petsc_level,grid[k].r);CHKERRQ(ierr);
    ierr = PCMGSetResidual(pc,petsc_level,PCMGDefaultResidual,grid[k].A); CHKERRQ(ierr);
//     ....Create interpolation between the levels....
    ierr = PCMGSetInterpolation(pc,petsc_level,grid[k].P);CHKERRQ(ierr);
    ierr = PCMGSetRestriction(pc,petsc_level,grid[k].R);CHKERRQ(ierr);  
// }
   /*....Set coarse grid solver....*/  
    ierr = PCMGGetCoarseSolve(pc,&coarsegridksp); CHKERRQ(ierr); 
    ierr = KSPSetFromOptions(coarsegridksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(coarsegridksp,grid[k+1].A,grid[k+1].A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(coarsegridksp,&coarsegridpc); CHKERRQ(ierr);
    
     k = k+1;
     if(k<=coarselevel){PetscPrintf(PETSC_COMM_WORLD,"\n  ------if works; coarsesolve with gmres for k = %d\n",k); 
//     ====================================================
//     Calls coarse solver recursively 
//     ====================================================
     ierr = KSPSetType(coarsegridksp, KSPFGMRES); CHKERRQ(ierr);
//      ierr = KSPGMRESSetRestart(coarsegridksp,grid[k].iter;) CHKERRQ(ierr); 
     ierr = KSPSetTolerances(coarsegridksp, 1e-12, 1e-50, 1e7,grid[k].iter);CHKERRQ(ierr); 
     ierr = PCSetType(coarsegridpc,PCSHELL); CHKERRQ(ierr); 
     ierr = MGSLPShellPCCreate(&shell); CHKERRQ(ierr);
     ierr = PCShellSetApply(coarsegridpc,MGSLPShellPCApply); CHKERRQ(ierr);
     ierr = PCShellSetContext(coarsegridpc,shell); CHKERRQ(ierr);
     ierr = PCShellSetDestroy(coarsegridpc,MGSLPShellPCDestroy);CHKERRQ(ierr);
     ierr = MGSLPShellPCSetUp(coarsegridpc,&grid[0],xforC,k); CHKERRQ(ierr);
     }
     else {PetscPrintf(PETSC_COMM_WORLD,"\n  ------else if; PCLU on coarsesolve for k = %d\n",k); 
     ierr = KSPSetType(coarsegridksp, KSPPREONLY); CHKERRQ(ierr);
     ierr = KSPSetTolerances(coarsegridksp, 1e-12, 1e-50, 1e7,1);CHKERRQ(ierr); 
     ierr = PCSetType(coarsegridpc, PCLU); CHKERRQ(ierr); 
     }
//       ierr = KSPView(coarsegridksp, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
   //   ierr = KSPSetMonitor(coarsegridksp,KSPDefaultMonitor,PETSC_NULL, 
   //          PETSC_NULL); CHKERRQ(ierr); 
 /*..Allow above criterea to be overwritten..*/
    ierr = PCSetFromOptions(pc); CHKERRQ(ierr);

// ierr = PCSetType(pc, PCILU); CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_SELF, "Inside user defined preconditioner\n");
    ierr = PCSetUp(pc); CHKERRQ(ierr);

    /*ierr = VecDestroy(&grid[k].b); CHKERRQ(ierr);
    ierr = VecDestroy(&grid[k].x); CHKERRQ(ierr);
    ierr = VecDestroy(&grid[k].r); CHKERRQ(ierr);
    ierr = VecDestroy(&grid[k].xforM); CHKERRQ(ierr);
    ierr = MatDestroy(&grid[k].A); CHKERRQ(ierr); 
    ierr = MatDestroy(&grid[k].M); CHKERRQ(ierr); 
    ierr = MatDestroy(&grid[k].P); CHKERRQ(ierr); 
    ierr = MatDestroy(&grid[k].R); CHKERRQ(ierr);*/ 
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
  PetscErrorCode ierr;
MGSLPShellPC *shell;

ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);
  /*..Free PCShell context..*/
  ierr = PetscFree(shell);CHKERRQ(ierr);

  return 0;
}