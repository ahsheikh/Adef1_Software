    
    ierr = VecCreate(PETSC_COMM_SELF,&tempvec); CHKERRQ(ierr);
    ierr = VecSetFromOptions(tempvec); CHKERRQ(ierr); 
    ierr = VecLoad(tempvec,fd);CHKERRQ(ierr);                      //  loads vector of same size as of no. of levels
    ierr = VecGetSize(tempvec,&lev); CHKERRQ(ierr);
    grid[0].lev = lev-1; 
    
//     ierr = VecView(tempvec,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);  // view vector
    
    ierr = MatCreate(PETSC_COMM_SELF,&(grid[1].A));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[1].A);CHKERRQ(ierr);
    ierr = MatLoad(grid[1].A,fd); CHKERRQ(ierr);            // loads fine grid helmholtz matrix. 
    
    ierr = MatCreate(PETSC_COMM_SELF,&(grid[1].M));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[1].M);CHKERRQ(ierr);
    ierr = MatLoad(grid[1].M,fd); CHKERRQ(ierr);          // loads find grid slp matrix
    
    
//     Vec ShiftVec; 
//     ierr = VecCreate(PETSC_COMM_SELF,&ShiftVec); CHKERRQ(ierr);
//     ierr = VecSetFromOptions(ShiftVec); CHKERRQ(ierr); 
//     ierr = VecLoad(ShiftVec,fd);CHKERRQ(ierr);
    

    for(i=1;i<lev+1;i++){
    //     ==========================================
    ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].P));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].P);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].P,fd); CHKERRQ(ierr);
    //     ==========================================
    ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].R));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].R);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].R,fd); CHKERRQ(ierr);
    //     ===========================================
    }
    
    ierr = VecCreate(PETSC_COMM_SELF,&b); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
    ierr = VecLoad(b,fd);CHKERRQ(ierr);
    
    
    // ========= creating A and M at levels other than fineest. 
    for(i=2;i<lev+2;i++){
      ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].A));CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].M));CHKERRQ(ierr);
    }
    
       
    for(i=1;i<lev+1;i++){                               // 2 ==> lev+1
      ierr = MatMatMult(grid[i].A, grid[i].P, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].A); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].R, grid[i+1].A, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].A); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].M, grid[i].P, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].M); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].R, grid[i+1].M, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].M); CHKERRQ(ierr);  
    }

     ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
     