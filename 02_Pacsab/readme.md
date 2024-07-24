# PACSAB
Molecular dynamics simulations of proteins are usually performed
on a single molecule, and coarse-grained protein models are calibrated
using single-molecule simulations. So there is a need for a model for the study of many protein systems, able to reproduce the properties of folded and unfolded proteins (in isolation and aggregation)
Standing for Pairwise Additive Potential for Coarse-Grained Side Chains and Atomistic Backbone, PACSAB is a protein model that simulates many-protein systems. It differs to other models by being based on contraction of a contraction of an implicit solvent classical atomistic model.

## Outputs from PACSAB
The simulation yields multiple snapshots, where each instance indicates the current state of the molecule. Here is an example of the snapshots produced by PACSAB:

**Snapshot 0**
|Atom string|Atom number|Atom name|Resiude name|Chain|Residue number|X coord|Y coord|Z coord|
|---|---|---|---|---|---|---|---|---|
|ATOM|1|N|ASP|A|1|153.059|135.000|152.027|
|ATOM|2|H|ASP|A|1|153.128|135.972|152.291|
|ATOM|3|CA|ASP|A|1|152.868|134.206|153.252|
|...|...|...|...|...|...|...|...|...|
|ATOM|1137|O|SER|B|89|154.974|154.441|125.225|
|ATOM|1138|S1|SER|B|89|155.855|157.643|124.970|
**Snapshot 261**
|Atom string|Atom number|Atom name|Resiude name|Chain|Residue number|X coord|Y coord|Z coord|
|---|---|---|---|---|---|---|---|---|
|ATOM|1|N|ASP|A|1|109.924|151.488|216.778|
|ATOM|2|H|ASP|A|1|106.467|151.607|216.072|
|ATOM|3|CA|ASP|A|1|110.739|150.099|216.317|
|...|...|...|...|...|...|...|...|...|
|ATOM|1137|O|SER|B|89|192.297|157.424|117.770|
|ATOM|1138|S1|SER|B|89|191.712|159.001|114.667|

Lastly, here is a visualization from these snapshots in a visual software (VMD 1.9.3)
![Simulation of two proteins in VMD](./simulation.gif)`
