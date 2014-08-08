ierr = VecCreate(PETSC_COMM_WORLD,&tempvec); CHKERRQ(ierr);
    ierr = VecSetFromOptions(tempvec); CHKERRQ(ierr); 
    ierr = VecLoad(tempvec,fd);CHKERRQ(ierr);
    ierr = VecGetSize(tempvec,&lev); CHKERRQ(ierr);
    grid[0].lev = lev;
     PetscPrintf(PETSC_COMM_WORLD,"\n  Total # of levels (at start)   =  %d\n",lev);
      
     for(i=1;i<lev+2;i++){
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].A));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].A);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].A,fd); CHKERRQ(ierr);
    
    //     ==========================================
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].M));CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].M));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].M);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].M,fd); CHKERRQ(ierr);
    //     ==========================================
//     ierr = MatView(grid[i].A,PETSC_NULL); CHKERRQ(ierr);
     }
//      for(i=1;i<5;i++){
//     ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].M));CHKERRQ(ierr);
//     ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].M));CHKERRQ(ierr);
//     ierr=  MatSetFromOptions(grid[i].M);CHKERRQ(ierr);
//     ierr = MatLoad(grid[i].M,fd); CHKERRQ(ierr);
// //     ierr = MatView(grid[i].M,PETSC_NULL); CHKERRQ(ierr);
//      }
     for(i=1;i<lev+1;i++){
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].P));CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].P));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].P);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].P,fd); CHKERRQ(ierr);
    
//     ==========================================
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].R));CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].R));CHKERRQ(ierr);
    ierr=  MatSetFromOptions(grid[i].R);CHKERRQ(ierr);
    ierr = MatLoad(grid[i].R,fd); CHKERRQ(ierr);
//     ===========================================
//     ierr = MatView(grid[i].P,PETSC_NULL); CHKERRQ(ierr);
     }
//      for(i=1;i<4;i++){
//     ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].R));CHKERRQ(ierr);
//     ierr = MatCreate(PETSC_COMM_WORLD,&(grid[i].R));CHKERRQ(ierr);
//     ierr=  MatSetFromOptions(grid[i].R);CHKERRQ(ierr);
//     ierr = MatLoad(grid[i].R,fd); CHKERRQ(ierr);
// //     ierr = MatView(grid[i].R,PETSC_NULL); CHKERRQ(ierr);
//      }
     
    ierr = MatGetSize(grid[2].A,PETSC_NULL,&grid[1].nH); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,&b); CHKERRQ(ierr);
    ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
    ierr = VecLoad(b,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);