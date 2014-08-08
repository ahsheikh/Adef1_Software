static char help[] = "Reads a PETSc 4 matrices A, AH, P, R, and a vector f \n\
and vector from a file and solves linear system Ax = f by Multigrid \n\
using KPS-richardson and MG as preconditioner. Thus cycles of Multigrid \n\
can be controlled by Richardson iteratiosn.  Input parameters include \n\
  -f <input_file> : file to load matrices and one vector \n\n";
// #include <math.h>
// #include <../src/mat/impls/aij/seq/aij.h>
// #include <../src/mat/impls/aij/mpi/mpiaij.h>
// #include <private/pcimpl.h>   /*I "petscpc.h" I*/
// #include <../src/ksp/pc/impls/mg/mgimpl.h>   

#include "petscsys.h" 
#include "petscksp.h"
#include "petscpcmg.h" 
#include "MGSLP.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
 
//   Mat            A,AH,P,R,FineLevelMatrix,M;            
  Vec            x,b,tempvec, xPC;          
  PC		 pc; 
  PetscViewer    fd;               /* viewer */
  char           file[PETSC_MAX_PATH_LEN];     /* input file name */
  KSP            ksp;        
  PetscBool      flg; 
  PetscErrorCode ierr;
  PetscMPIInt    lev, rank=4;
  PetscLogDouble v1,v2,t_setup, t_solve;
  PetscInt 	 its,i, ii,k;
  MGSLPShellPC 	 *shell;
//   MGforM         *shellforM;
  GridCtx	 grid[MAX_LEVELS]; 
  
  
  /*..Executable statements..*/
    PetscInitialize(&argc,&args,(char *)0,help);
  
  /*..Get start time of linear system setup..*/ 
    ierr = PetscGetTime(&v1);CHKERRQ(ierr); 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  /*...  Determine files from which we read the operators ... */
    PetscOptionsGetString(PETSC_NULL,"-f",file,PETSC_MAX_PATH_LEN-1,&flg);
  /* Open binary file.  Note that we use FILE_MODE_READ to indicate reading from this file. */
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd); CHKERRQ(ierr);
  /* Creat and  Load the matrices and vector; then destroy the viewer.  */
// #include "dataload.h"
// #include "dataloadADV.h"
#include "loadingData.h"
int tempsizeA; 
  ierr = MatGetSize(grid[1].A, &tempsizeA, &tempsizeA); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\n  ------size of fine matrix A = %d\n",tempsizeA); 
//      Creating solution vector adn naming it. 
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr); ierr = VecDuplicate(b,&xPC);CHKERRQ(ierr);
//      CREATING KSP SOLVER CONTEXT 
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
//      Setting operators. Here the matrix that defines the linear system also serves as the preconditioning matrix.
    ierr = KSPSetOperators(ksp,grid[1].A,grid[1].M,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
//      KSP is specified with some options. These can be overruled by 
//      giving different in runtime. This is faciliates by KSPSetFromOptions      
     ierr = KSPSetTolerances(ksp,1.e-6,1.e-50,1.e5,PETSC_DEFAULT);CHKERRQ(ierr);
     ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
     ierr = KSPSetType(ksp,KSPFGMRES); CHKERRQ(ierr); 
//      ierr = KSPGCRSetRestart(ksp,100);CHKERRQ(ierr);
     ierr = KSPGMRESSetRestart(ksp,500); CHKERRQ(ierr); 
//      ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 

//     Set runtime options, e.g.,
//         -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
//     These options will override those specified above as long as
//     KSPSetFromOptions() is called _after_ any other customization routines.
     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
     
//   =======================================================================
//     Here we call user-defined preconditioner  
//      ierr = PCSetType(pc,PCLU); CHKERRQ(ierr); 
//   =======================================================================

     ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr); 
     ierr = MGSLPShellPCCreate(&shell); CHKERRQ(ierr);
     ierr = PCShellSetApply(pc,MGSLPShellPCApply); CHKERRQ(ierr);
     ierr = PCShellSetContext(pc,shell); CHKERRQ(ierr);
     ierr = PCShellSetDestroy(pc,MGSLPShellPCDestroy);CHKERRQ(ierr);
     ierr = PCShellSetName(pc,"MyPreconditioner");CHKERRQ(ierr); 
     
          k = 1; 
     grid[2].iter = 8; grid[3].iter = 4; grid[4].iter = 1; grid[5].iter = 1;
     grid[6].iter = 1; grid[7].iter = 1; grid[7].iter = 1; grid[8].iter = 1;
     grid[9].iter = 1; grid[10].iter = 1; grid[11].iter = 1;
     
     
     int IterLev5, IterLev6, IterLev7, IterLev8, IterLev9; 
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev2", &grid[2].iter, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev3", &grid[3].iter, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev4", &grid[4].iter, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev5", &IterLev5, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev6", &IterLev6, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev7", &IterLev7, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev8", &IterLev8, 0);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL, "-iter_lev9", &IterLev9, 0);CHKERRQ(ierr);    

     ierr = MGSLPShellPCSetUp(pc,&grid[0],xPC,k); CHKERRQ(ierr);
         
//   ======================================================================

     
       ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
//      Indicate that we are going to use a non-zero initial solution
//      ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);

//      Get end time of linear system setup
     ierr = PetscGetTime(&v2);CHKERRQ(ierr); 
     t_setup = v2 - v1;  

//      Get start time of linear solve
    ierr = PetscGetTime(&v1);CHKERRQ(ierr); 
 
//      ierr = SLESView(sles, VIEWER_STDOUT_SELF); CHKERRQ(ierr);
//      Solve linear system
    ierr = KSPSolve(ksp,b,x); CHKERRQ(ierr);
//      ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
         ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);

//      Print solution vector
//      ierr = VecView(x,PETSC_NULL); CHKERRQ(ierr);
    
//      Print number of iterations
    PetscPrintf(PETSC_COMM_WORLD,"\n** Number of iterations done = %d \n",its);   

  /*..Get end time of linear solve..*/ 
    ierr = PetscGetTime(&v2);CHKERRQ(ierr); 
    t_solve = v2 - v1;  

   printf("\n[PETSc]:Time spend in setup = %e \n",t_setup);
   printf("[PETSc]:Time spend in solve = %e \n",t_solve);
   printf("[PETSc]:Total time = %e \n\n", t_setup + t_solve);
   
   ierr = VecDestroy(&x); CHKERRQ(ierr);
   ierr = VecDestroy(&b); CHKERRQ(ierr); 
//    ierr = MatDestroy(&grid[1].A); CHKERRQ(ierr); 
   ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
//    ierr = VecDestroy(&tempvec); CHKERRQ(ierr);
   ierr = VecDestroy(&xPC); CHKERRQ(ierr);
   
   for(ii=1;ii<lev+1;ii++){
    ierr = MatDestroy(&(grid[ii].P));CHKERRQ(ierr);
    ierr = MatDestroy(&(grid[ii].R));CHKERRQ(ierr);
    ierr = MatDestroy(&(grid[ii].A));CHKERRQ(ierr);
    ierr = MatDestroy(&(grid[ii].M));CHKERRQ(ierr);
  }
//    ierr = MatDestroy(&(grid[lev+2].A));CHKERRQ(ierr);
//    ierr = MatDestroy(&(grid[lev+2].M));CHKERRQ(ierr); 
   
   /*  Next function allows to give an other pc type at runtime */
  //   ierr = PCSetFromOptions(pc);CHKERRQ(ierr);
       /* Solve linear system  */
     /* ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);*/
    
   /* FINALIZING PETSC CODE 	*/
     ierr = PetscFinalize();CHKERRQ(ierr);
     return 0;  
 }
