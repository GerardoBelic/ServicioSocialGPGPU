# Refactor of PACSAB code
Just one step before translating the C++ code into Julia, I'm going to separate some code into functions and delete some unnecessary nested code.
## List of changes
- Removal of constants that limited the molecules maximum size, effectively making the program capable of accepting any group of molecules.
- Commented all unused variables declared in the structs used in the program (Xoc, Atpres, Intr, ...).
- Simplification of nested code, making the code easier to read.
- Separation of the initial part of the code intro functions. Each function has a specific purpose like fetching the molecules info, get the atom types or parameters, etc. I left the loop part of the code untouched because its in our interest making this part of the code faster (the first part of the program is fast but is not important to make it faster).

## Compilation and execution

    cd /02_Refactor
    
    # C++ compilation
    g++ -I"../../dependencies/include" Cpp/cgdmd.cpp -o Cpp/cgdmd
    
    # To test C++ (cut execution with Ctrl+C at any time)
    .\Cpp\cgdmd
    
This should produce the exact same output as the previous step (01_Translate_&_Debugging)
