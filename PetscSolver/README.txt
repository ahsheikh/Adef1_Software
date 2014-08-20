
Note that main program MainSolver.c make call of two other functions;

    1. MLadef.c
    2. GMGcycle.c

Write these files accordingly in your make file and 
keep that make file in directory, where main program MainSolver.c is . 


==============================


An example of how to compile and execute the main file and sub files. 

> make GMGcycle.c; make MLadef.c; make MainSolver

./MainSolver -f ~/path/to/your/datafile/filename.dat 