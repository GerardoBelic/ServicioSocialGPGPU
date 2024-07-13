# Composition of PACSAB
PACSAB is made up of 3 programs:

 - ***gentop.f***: takes a molecule as its input and groups some elements into one single "element". The output is another molecule with less size (due to the elements grouping) and a topology file with the distances between elements (the ones greater than 1.0 angstrom).
 - ***gentopcg.f***: takes the molecule and the topology files produced by *gentop.f* and groups some smaller groups into bigger ones. The output is again a molecule with less size and a topology file with the distances of the elements/groups that meet some criteria.
 - ***cgdmd.f***: the main file, takes the molecule and the topology files produced by *gentopcg.f* and make the protein/molecule unfolding. It produces many outputs, the main one being a snapshots file containing every step the molecule did during the unfolding.

## Compilation and execution

    cd /00_Pacsab_Fortran
    
    gfortran -ffree-line-length-0 pacsab/gentop.f95 -o pacsab/gentop
    gfortran -ffree-line-length-0 pacsab/gentopcg.f95 -o pacsab/gentopcg
    gfortran -ffree-line-length-0 pacsab/cgdmd.f95 -o pacsab/cgdmd
    
    # Pick a molecule from /Test_Molecules, for example leapoutstruct1AY7
    .\"pacsab/gentop" < ../Test_Molecules/leapoutstruct1AY7.pdb
    .\"pacsab/gentopcg" < structure.pdb
    .\"pacsab/cgdmd" < input/dmdcg.dat
    
    # The execution takes some time, so stop the execution at any time with Ctrl+C

The main output of this file is snapshots.pdb, which holds every iteration of the unfolding. We can visualize graphically every frame with VMD or a similar software.

