 
    ierr = VecCreate(PETSC_COMM_SELF,&tempvec); CHKERRQ(ierr);
    ierr = VecSetFromOptions(tempvec); CHKERRQ(ierr); 
    ierr = VecLoad(tempvec,fd);CHKERRQ(ierr);
    ierr = VecGetSize(tempvec,&lev); CHKERRQ(ierr);
    grid[0].lev = lev;
    
    ierr = VecView(tempvec,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
    ierr = MatCreate(PETSC_COMM_SELF,&(grid[1].A));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[1].A);CHKERRQ(ierr);
    ierr = MatLoad(grid[1].A,fd); CHKERRQ(ierr);
    Vec ShiftVec; 
    ierr = VecCreate(PETSC_COMM_SELF,&ShiftVec); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ShiftVec); CHKERRQ(ierr); 
    ierr = VecLoad(ShiftVec,fd);CHKERRQ(ierr);
    

    for(i=1;i<lev+1;i++){
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
    
    

// =============== SETTING UP A & M matrices ============================    
    PetscScalar    shift_a, shift_m, zero = 0;
    PetscInt 	   sizeVec; 
    Vec 	   yVec, da, dm; 
    PetscRealPart(shift_a)=1; PetscImaginaryPart(shift_a)=0;
    PetscRealPart(shift_m)=1; PetscImaginaryPart(shift_m)=-1;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-shift_a",&shift_a,PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-shift_m",&shift_m,PETSC_NULL); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_SELF, &yVec); CHKERRQ(ierr);
    ierr = VecSetFromOptions(yVec);CHKERRQ(ierr);
    ierr = VecGetSize(ShiftVec,&sizeVec); CHKERRQ(ierr); 
    ierr = VecSetSizes(yVec,PETSC_DECIDE,sizeVec); CHKERRQ(ierr); 
    ierr = VecSet(yVec,zero);CHKERRQ(ierr);
    
    ierr = VecDuplicate(yVec,&da); CHKERRQ(ierr);
    ierr = VecDuplicate(yVec,&dm); CHKERRQ(ierr);     

    ierr = VecWAXPY(da, shift_a,ShiftVec,yVec); CHKERRQ(ierr);
    ierr = VecWAXPY(dm, shift_m,ShiftVec,yVec); CHKERRQ(ierr); 
    
    ierr = MatDuplicate(grid[1].A,MAT_COPY_VALUES,&grid[1].M); CHKERRQ(ierr);
    ierr = MatCopy(grid[1].A,grid[1].M,SAME_NONZERO_PATTERN);CHKERRQ(ierr); 

    for(i=2;i<lev+2;i++){
      ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].A));CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_SELF,&(grid[i].M));CHKERRQ(ierr);
    }
    
    ierr = MatDiagonalSet(grid[1].A,da,ADD_VALUES); CHKERRQ(ierr);
    ierr = MatDiagonalSet(grid[1].M,dm,ADD_VALUES); CHKERRQ(ierr); 
    
    for(i=1;i<lev+1;i++){                               // 2==> lev+1
      ierr = MatMatMult(grid[i].A, grid[i].P, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].A); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].R, grid[i+1].A, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].A); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].M, grid[i].P, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].M); CHKERRQ(ierr);
      ierr = MatMatMult(grid[i].R, grid[i+1].M, MAT_INITIAL_MATRIX, 1.0, &grid[i+1].M); CHKERRQ(ierr);  
    }

     ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
     ierr = VecDestroy(&da); CHKERRQ(ierr);
     ierr = VecDestroy(&dm); CHKERRQ(ierr);
     ierr = VecDestroy(&yVec); CHKERRQ(ierr);
     ierr = VecDestroy(&ShiftVec); CHKERRQ(ierr);
     ierr = VecDestroy(&tempvec); CHKERRQ(ierr);