# Gro: set of classes to handle Gromacs input and output data

## Classes

### GMolecule

Represents an individual molecule. 
- **Z**: atomic numbers;
- **xyz**: coordinates in Å;
- **name**: short name of the compound. Typically, for Gromacs studied species is `LIG`, and solvent is `SOL` (water can be `WAT` as well);
- **index**: unique index of the molecule.

### GFrame

Represents an individual snapshot from MD simulation, i.e., a collection of *GMolecule*s. Can handle a single \*.gro file (which can be an input for a simulation run as well as an output).

### GTrajectory

Represents a MD trajectory, i.e. a collection of *GFrame*s. Can handle raw binary \*.trr file and also ASCII \*.gro file.

## Features

- Automatically fixes the molecules at the edge of the simulation box. Typically raw output from Gromacs  contains 'broken' molecules, e.g. some atoms are on another side of the box than the others. Currently `fix()` is called for each molecule which makes reading of large TRR files rather slow. To be optimized.
- Deduces the local neighbourhood for selected molecules. This is used to generate an input for e.g. spectral calculations, which require only a small blob of atoms. 
- Removes virtual atoms, used by  some MD simulations E.g. widely used TIP4P model of water consists of 4 atoms (O, H, H and one virtual). 