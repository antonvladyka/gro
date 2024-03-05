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

## Examples

Individual gromacs snapshot

    # read individual snapshot
    filename = 'md.gro'
    gro = GFrame(filename)

    # molecule in the frame
    print(gro.n_molecules_)

    # xyz array
    print(gro.xyz)

    # get neighbourhood for a given molecule
    ## the following code returns the system of all the molecules, whoch are neighbours of molecule #10
    ## neighbours are the molecules such that all atoms are within given distance threshold from the target
    ## resulting system is centered at atom #6 of the molecule #10 (convenient for visualization or spectral
    ## calculations)
    nei = gro.neighbourshood(10, dist={'SOL': 5, 'LIG': 7.5}, center_at=6)
    nei.to_xyz() # prints the content of the system
    nei.to_xyz(filename=<>) # saves .xyz file

Trajectory (requires also gro file, since trajectorry file containes only coordinates, timing, box size, but ton atomic numbers)
    
    file_trr = 'md.trr'
    file_gro = 'md.gro'
    traj = GTrajectory(file_trr, file_gro)  # takes a while to read binary data and also fix the molecules
    gro = traj[0]  # get the first frame of the trajectory
