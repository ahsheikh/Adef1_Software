   
    int k,petsc_level,  levels,levelsforM,sizeA;
    KSP coarsegridksp; 
    PC  coarsegridpc; 

    ierr = PCSetType(pc,PCMG);CHKERRQ(ierr);
    levels = grid[0].lev+1;
    ierr = PetscOptionsGetInt(PETSC_NULL, "-levels", &levels, 0);CHKERRQ(ierr);
    ierr = PCMGSetLevels(pc,levels,PETSC_NULL);CHKERRQ(ierr);

for (k=1; k<levels; k++){  
     petsc_level = levels - k;  /*petsc_level is current level.dont give coarsest = 0.*/
     levelsforM = petsc_level+1; 
//     ....Get pre-smoothing context....
    /*ierr = PCMGGetSmootherDown(pc,petsc_level,&grid[k].ksp_pre); CHKERRQ(ierr); 
//     ierr = KSPView(grid[k].ksp_pre,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
    ierr = PCMGGetSmootherUp(pc,petsc_level,&grid[k].ksp_post); CHKERRQ(ierr); 
  */ 
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
//     ierr = VecDuplicate(grid[k].x,&grid[k].b); CHKERRQ(ierr); 
//     ierr = VecDuplicate(grid[k].x,&grid[k].r); CHKERRQ(ierr);

//     ...set ksp_pre context....
//     ierr = KSPSetOperators(grid[k].ksp_pre,grid[k].A,grid[k].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//     ierr = KSPGetPC(grid[k].ksp_pre,&grid[k].pc_pre);CHKERRQ(ierr);
//     ierr = KSPSetType(grid[k].ksp_pre, KSPPREONLY);CHKERRQ(ierr);
// //     ierr = KSPRichardsonSetScale(grid[k].ksp_pre, omega); CHKERRQ(ierr);
//     ierr = KSPSetTolerances(grid[k].ksp_pre, 1e-12, 1e-50, 1e7,0);CHKERRQ(ierr); 
//     PetscPrintf(PETSC_COMM_SELF, "Inside user defined preconditioner before applying MG on M \n"); 
//     ierr = PCSetType(grid[k].pc_pre, PCNONE);CHKERRQ(ierr); 
// ================================ TO CHECK an other PCSHELL_FOR_M which applies MG on M instead LU =======
//      ierr = PCSetType(grid[k].pc_pre,PCSHELL); CHKERRQ(ierr); 
//      ierr = MGforMCreate(&shellforM); CHKERRQ(ierr);
//      ierr = PCShellSetApply(grid[k].pc_pre,MGforMApply); CHKERRQ(ierr);
//      ierr = PCShellSetContext(grid[k].pc_pre,shellforM); CHKERRQ(ierr);
//      ierr = PCShellSetDestroy(grid[k].pc_pre,MGforMDestroy);CHKERRQ(ierr);
//      ierr = PCShellSetName(grid[k].pc_pre,"MyPreconditionerforM");CHKERRQ(ierr); 
//      ierr = VecDuplicate(grid[k].b,&xforM);CHKERRQ(ierr);
//      ierr = MGforMSetUp(grid[k].pc_pre,&grid[0],xforM,levelsforM); CHKERRQ(ierr);
// =========================================================================================================
//     ....set KSP_POST context....
//     ierr = KSPSetOperators(grid[k].ksp_post,grid[k].A,grid[k].A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//     ierr = KSPGetPC(grid[k].ksp_post,&grid[k].pc_post); CHKERRQ(ierr);
//     ierr = KSPSetInitialGuessNonzero(grid[k].ksp_post, PETSC_TRUE);CHKERRQ(ierr);
//     ierr = KSPSetType(grid[k].ksp_post,KSPRICHARDSON); CHKERRQ(ierr);
//     ierr = KSPRichardsonSetScale(grid[k].ksp_post, omega); 
//     ierr = KSPSetTolerances(grid[k].ksp_post,1e-12, 1e-50, 1e7,0); CHKERRQ(ierr);  
//     ierr = PCSetType(grid[k].pc_post, PCNONE); CHKERRQ(ierr); 
//     ierr = PCSetType(grid[k].pc_post, PCKSP); CHKERRQ(ierr); 
//     ierr = PCKSPGetKSP(grid[k].pc_post,&ksp1); CHKERRQ(ierr);
//     ierr = KSPSetType(ksp1, KSPPREONLY);CHKERRQ(ierr);
//     ierr = KSPGetPC(ksp1, &pc1);CHKERRQ(ierr);
//     ierr = PCSetType(pc1,PCNONE);CHKERRQ(ierr);

// ================================ TO CHECK an other PCSHELL_FOR_MULTIPLICATION_OF_A  ================
//      ierr = PCSetType(grid[k].pc_post,PCSHELL); CHKERRQ(ierr); 
//      ierr = pcshiftCreate(&shellpcshift); CHKERRQ(ierr);
//      ierr = PCShellSetApply(grid[k].pc_post,pcshiftApply); CHKERRQ(ierr);
//      ierr = PCShellSetContext(grid[k].pc_post,shellpcshift); CHKERRQ(ierr);
//      ierr = PCShellSetDestroy(grid[k].pc_post,pcshiftDestroy);CHKERRQ(ierr);
//      ierr = PCShellSetName(grid[k].pc_post,"MyPreconditionerforMultiplicationOfA+Shift");CHKERRQ(ierr); 
//      ierr = VecDuplicate(grid[k].b,&x_pcshift);CHKERRQ(ierr);
//      ierr = pcshiftSetUp(grid[k].pc_post,&grid[0],x_pcshift); CHKERRQ(ierr);
// =========================================================================================================
//     ============================================
    
    ierr = PCMGSetX(pc,petsc_level-1,grid[k].x);CHKERRQ(ierr); 
    ierr = PCMGSetRhs(pc,petsc_level-1,grid[k].b);CHKERRQ(ierr); 
//     ierr = PCMGSetR(pc,petsc_level,grid[k].r);CHKERRQ(ierr);
    ierr = PCMGSetResidual(pc,petsc_level,PCMGDefaultResidual,grid[k].A); CHKERRQ(ierr);
//     ....Create interpolation between the levels....
    ierr = PCMGSetInterpolation(pc,petsc_level,grid[k].P);CHKERRQ(ierr);
    ierr = PCMGSetRestriction(pc,petsc_level,grid[k].R);CHKERRQ(ierr);  
//      ierr = PCMGSetInterpolation(pc,1,R);CHKERRQ(ierr);
//      ierr = PCMGSetRestriction(pc,1,P);CHKERRQ(ierr);  
}
  /*....Set coarse grid solver....*/  
    ierr = PCMGGetCoarseSolve(pc,&coarsegridksp); CHKERRQ(ierr); 
    ierr = KSPSetFromOptions(coarsegridksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(coarsegridksp,grid[levels].A,grid[levels].A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    ierr = KSPGetPC(coarsegridksp,&coarsegridpc); CHKERRQ(ierr);
    ierr = KSPSetType(coarsegridksp, KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPSetTolerances(coarsegridksp, 1e-12, 1e-50, 1e7,1);CHKERRQ(ierr); 
    ierr = PCSetType(coarsegridpc, PCLU); CHKERRQ(ierr); 
   //   ierr = KSPView(coarsegridksp, VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
   //   ierr = KSPSetMonitor(coarsegridksp,KSPDefaultMonitor,PETSC_NULL, 
   //          PETSC_NULL); CHKERRQ(ierr); 
 /*..Allow above criterea to be overwritten..*/
    ierr = PCSetFromOptions(pc); CHKERRQ(ierr);