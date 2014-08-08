//  Here we define user-defined function which applied Multigrid on M
#include "MGSLP.h"
#include "petscksp.h"
#include "petscpcmg.h"
// ===================================================================================================================
//                         here we start to defined MG approximation for M 
// ===================================================================================================================

/* ------------------------------------------------------------------- */
#undef __FUNC__
#define __FUNC__ "MGSforMCreate"

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
 PetscErrorCode MGforMSetUp(PC pcforM,GridCtx* grid,Vec x,PetscInt levelforM)
 {
    PetscErrorCode ierr;
    int ind, petsc_levelforM, levels, sizeM, levelsforM;
    PetscReal    omega;
    KSP  CGkspForM; 
    PC   CGpcForM; 
    omega = 0.7;               /*Relaxation factor*/
PetscPrintf(PETSC_COMM_WORLD,"Level for M = %d\n",levelforM);
//     ierr = PCSetType(pcforM,PCLU);CHKERRQ(ierr);
//     PetscPrintf(PETSC_COMM_SELF, "Here MGforM applies PCLU on M \n");
// ================================PCMG  structure for M ===============================================
//     PetscPrintf(PETSC_COMM_SELF, "Here MGforM prepares to apply Vcycle on M for k =  %d\n",levelforM);     

    ierr = PCSetType(pcforM,PCMG);CHKERRQ(ierr);
    levels = grid[0].lev+1;  
//     PetscPrintf(PETSC_COMM_SELF, "tells total number of levels inside MG approximation of M %d\n",levels);
    levelsforM = levels-levelforM +1;
    
    ierr = PetscOptionsGetInt(PETSC_NULL, "-levelsforM", &levelsforM,0);CHKERRQ(ierr);
//      PetscPrintf(PETSC_COMM_SELF, "Tells NEW GIVEN number of levels inside MGforM on M =  %d\n",levelsforM);
//      PetscPrintf(PETSC_COMM_SELF, "and starting point is k =  %d\n",levelforM);     
    
    ierr = PCMGSetLevels(pcforM,levelsforM,PETSC_NULL);CHKERRQ(ierr);
 
//     ierr = PCMGSetType(pcforM, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
//     ierr = PCMGMultiplicativeSetCycles(pcforM,1);CHKERRQ(ierr);
    
    ierr = PCMGSetType(pcforM, PC_MG_FULL); CHKERRQ(ierr);
//     ierr = PCMGSetNumberSmoothDown(pcforM,1); CHKERRQ(ierr);
//     ierr = PCMGSetNumberSmoothUp(pcforM,1); CHKERRQ(ierr);
for (ind=levelforM; ind<levels; ind++){  
     petsc_levelforM = levels - ind;  /*petsc_level is current level.dont give coarsest = 0.*/
//      PetscPrintf(PETSC_COMM_SELF, "Prints each level inside MGforM  %d\n",petsc_levelforM);     
//     ....Get pre-smoothing context....
    ierr = PCMGGetSmootherDown(pcforM,petsc_levelforM,&grid[ind].ksp_preForM); CHKERRQ(ierr); 
//     ierr = KSPView(grid[ind].ksp_preForM,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
    ierr = PCMGGetSmootherUp(pcforM,petsc_levelforM,&grid[ind].ksp_postForM); CHKERRQ(ierr); 
//    
    ierr = MatGetSize(grid[ind+1].M, &sizeM, &sizeM); CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].x_M); CHKERRQ(ierr);
    ierr = VecSetType(grid[ind].x_M, VECSEQ); CHKERRQ(ierr);
    ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].b_M); CHKERRQ(ierr);
    ierr = VecSetType(grid[ind].b_M, VECSEQ); CHKERRQ(ierr); 
    ierr = VecCreate(MPI_COMM_WORLD,&grid[ind].r_M); CHKERRQ(ierr);
    ierr = VecSetType(grid[ind].r_M,VECSEQ); CHKERRQ(ierr);
    ierr = VecSetSizes(grid[ind].x_M,PETSC_DECIDE,sizeM); CHKERRQ(ierr);
    ierr = VecSetSizes(grid[ind].r_M,PETSC_DECIDE,sizeM);  CHKERRQ(ierr);
    ierr = VecSetSizes(grid[ind].b_M,PETSC_DECIDE,sizeM);  CHKERRQ(ierr);
//     ierr = VecDuplicate(grid[k].x,&grid[k].b); CHKERRQ(ierr); 
//     ierr = VecDuplicate(grid[k].x,&grid[k].r); CHKERRQ(ierr);
//     PetscPrintf(PETSC_COMM_SELF, "size of M to be used as next level vector's size = %d\n",sizeM); 

//     ...set ksp_preForM context....
    ierr = KSPSetOperators(grid[ind].ksp_preForM,grid[ind].A,grid[ind].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(grid[ind].ksp_preForM,&grid[ind].pc_preForM);CHKERRQ(ierr);
    ierr = KSPSetType(grid[ind].ksp_preForM, KSPRICHARDSON);CHKERRQ(ierr);
   //ierr = KSPSetTolerances(grid[ind].ksp_preForM, 1e-12, 1e-50, 1e7,1);CHKERRQ(ierr); 
//     PetscPrintf(PETSC_COMM_SELF, "Performing PREsmoothing with Damped Jacobi on M @ Level = %d\n",petsc_levelforM); 
    ierr = KSPRichardsonSetScale(grid[ind].ksp_preForM, omega); CHKERRQ(ierr);
    ierr = PCSetType(grid[ind].pc_preForM, PCJACOBI);CHKERRQ(ierr); 
    ierr = PCMGSetNumberSmoothDown(pcforM,1); CHKERRQ(ierr);
//     ....set KSP_POST context....
    ierr = KSPSetOperators(grid[ind].ksp_postForM,grid[ind].M,grid[ind].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(grid[ind].ksp_postForM,&grid[ind].pc_postForM); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(grid[ind].ksp_postForM, PETSC_TRUE);CHKERRQ(ierr);
    ierr = KSPSetType(grid[ind].ksp_postForM, KSPRICHARDSON); CHKERRQ(ierr);
  //  ierr = KSPSetTolerances(grid[ind].ksp_postForM,1e-12, 1e-50, 1e7,1); CHKERRQ(ierr);  
//     PetscPrintf(PETSC_COMM_SELF, "Performing POSTsmoothing with Damped Jacobi on M @ Level = %d\n",petsc_levelforM); 
    ierr = KSPRichardsonSetScale(grid[ind].ksp_postForM, omega); CHKERRQ(ierr);
    ierr = PCSetType(grid[ind].pc_postForM, PCJACOBI); CHKERRQ(ierr);  
    ierr = PCMGSetNumberSmoothUp(pcforM,1); CHKERRQ(ierr);
//     ============================================
    ierr = PCMGSetX(pcforM,petsc_levelforM-1,grid[ind].x_M);CHKERRQ(ierr); 
    ierr = PCMGSetRhs(pcforM,petsc_levelforM-1,grid[ind].b_M);CHKERRQ(ierr); 
//     ierr = PCMGSetR(pc,petsc_level,grid[k].r);CHKERRQ(ierr);
    ierr = PCMGSetResidual(pcforM,petsc_levelforM,PCMGDefaultResidual,grid[ind].M); CHKERRQ(ierr);
//     ....Create interpolation between the levels....
    ierr = PCMGSetInterpolation(pcforM,petsc_levelforM,grid[ind].P);CHKERRQ(ierr);
    ierr = PCMGSetRestriction(pcforM,petsc_levelforM,grid[ind].R);CHKERRQ(ierr);  
//      ierr = PCMGSetInterpolation(pc,1,R);CHKERRQ(ierr);
//      ierr = PCMGSetRestriction(pc,1,P);CHKERRQ(ierr);  
}
  /*....Set coarse grid solver....*/  
    ierr = PCMGGetCoarseSolve(pcforM,&CGkspForM); CHKERRQ(ierr); 
    ierr = KSPSetFromOptions(CGkspForM); CHKERRQ(ierr);
    ierr = KSPSetOperators(CGkspForM,grid[levels].M,grid[levels].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(CGkspForM,&CGpcForM); CHKERRQ(ierr);
    ierr = KSPSetType(CGkspForM, KSPPREONLY); CHKERRQ(ierr);
    ierr = PCSetType(CGpcForM, PCLU); CHKERRQ(ierr); 
//      ierr = KSPView(coarsegridksp, VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
//      ierr = KSPSetMonitor(coarsegridksp,KSPDefaultMonitor,PETSC_NULL, 
//             PETSC_NULL); CHKERRQ(ierr); 
 /*..Allow above criterea to be overwritten..*/
    ierr = PCSetFromOptions(pcforM); CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_SELF, "Inside user defined preconditioner\n");
// =====================================================================================================
    ierr = PCSetUp(pcforM); CHKERRQ(ierr);
    for (ind=levelforM; ind<levels; ind++){  
       ierr = VecDestroy(&grid[ind].b_M); CHKERRQ(ierr);
       ierr = VecDestroy(&grid[ind].x_M); CHKERRQ(ierr);
       ierr = VecDestroy(&grid[ind].r_M); CHKERRQ(ierr);
       
    }
//     ierr = MatDestroy(&grid[ind].R); CHKERRQ(ierr);
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
//  
#undef __FUNC__
#define __FUNC__ "MGforMDestroy"
 PetscErrorCode MGforMDestroy(PC pcforM) 
{
  PetscErrorCode ierr;
  /*..Free PCShell context..*/

//    PetscFree(pcforM);

MGforM *shellforM;

ierr = PCShellGetContext(pcforM,(void**)&shellforM); CHKERRQ(ierr);
// ierr = VecDestroy(&shell->diag); CHKERRQ(ierr);
ierr = PetscFree(shellforM); CHKERRQ(ierr); 
//    ierr = PCDestroy(&pcforM); CHKERRQ(ierr); 
return 0;
}
