# Translation into C/C++
Before we can parallelize the program into CUDA, we need to translate the program into a more modern language (C/C++ and Julia).
It is planned to make the parallelization in Julia, but for now I'm more familiar with C++ so I'm making the translation in it. Because we need to separate the program's logic into functions (in the FORTRAN code there is almost no functions and a lot of repeated code), it doesn't matter much which language we choose, the final program will not have many lines so the translation from C++ to Julia should be no problem; plus we can compare the performance between FORTRAN, C++ and Julia.

## Goals
Translate the FORTRAN main code *cgdmd.f95* to C++. The two programs need to have the EXACT same output so the C++'s version can be considered valid.
Some things to be taken in mind:

 - In one section of *cgdmd.f95* (lines 984-996), we select a random atom and randomize its velocity. So we make `iterm = 0` in the original and translated codes to control the stochastic nature of the program.
 - Some other parts of the program use random numbers between 0 and 1 to assign values. To control this values, a new function was created that returns fake random numbers (in both FORTRAN and in C++ versions).
 - To get the exact same values (and by exact I mean exactly the same), I had to modify the FORTRAN constants/literals to be of type `real*8`, this is because the assignation `rohmin = 1.75` is different to `rohmin = 1.75d0`, the first one is declared as a `real*4/float`, and the second one is a `real*8/double` (yielding slightly but different results). Another change was the `sqrt(x)` functions in FORTRAN, for some reason they perform different than the operation `x**0.5d0`, the last one being the equivalent to C++'s `std::sqrt(x)`.

## Considerations
The source files in this path have a lot of commented code due to the debugging, so they are not very readable in this stage.

## Compilation and execution

    cd /01_Translate_&_Debugging
    
    # FORTRAN compilation
    gfortran -ffree-line-length-0 FORTRAN/cgdmd.f95 -o FORTRAN/cgdmd
    
    # C++ compilation
    g++ -I"../../dependencies/include" Cpp/cgdmd.cpp -o Cpp/cgdmd
    
    # To test FORTRAN (cut the program's execution at any time with Ctrl+C)
    .\FORTRAN\cgdmd < dmdcg.dat
    
    # To test C++ (same, cut execution with Ctrl+C at any time)
    .\Cpp\cgdmd
    
   Both programs will produce the same output.
