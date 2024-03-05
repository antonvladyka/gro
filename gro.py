import numpy as np
import re
import struct
from scipy.spatial.distance import cdist
from dataclasses import dataclass

_ATOMS = {'C': 6, 'H': 1, 'N': 7, 'O' : 8, 'OW': 8, 'HW': 1, 'MW': 0}  # MW is a virtual atom from tip4p water model. excluded via Molecule.fix() method
_CHARGES = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
    
@dataclass
class GMolecule:
    xyz: np.array
    z: np.array
    name: str
    index: int
    
    @property
    def fullname(self):
        return f'{self.name}{self.index}'
    
    def __repr__(self):
        return f'{self.fullname} ({len(self.z)} atoms)'
    
    def fix(self, box):
        """ 1) Excludes virtual atoms with z == 0
            2) Merges the broken atoms together, e.g., if atoms of one molecules are on other sides of the simulation box"""
        sub = self.z != 0
        self.xyz = self.xyz[sub]
        self.z = self.z[sub]
        xyz0 = self.xyz[0]
        d = np.round((self.xyz - xyz0)/box)
        if d.any():
            self.xyz -= d * box
        
    def center(self, xyz0, box):
        """Moves the molecule in a postion if the box is in centers at xyz0 point"""
        xyz = self.xyz
        d = np.round((xyz[0] - xyz0)/box)
        return type(self)(xyz - xyz0 - d * box, self.z, self.name, self.index)
    

class GFrame:
    """
    _parse : file/file descriptor -> list atoms
    _build_molecules
    TODO: implement read_legacy method to read fixed-width GRO format (can be also faster)
    """
    def __init__(self, param=None):
        self.molecules = []
        self.indices = []
        self.box = None
        self.n_molecules = 0
        self.t = None
        if param is None:
            return
        if type(param) == str:
            self.from_file(param)
        if type(param) == list:
            self.from_molecules(param)
    
    def __getitem__(self, idx):
        return self.molecules[idx]
    
    @classmethod
    def _parse(cls, f):
        """Reads opened *.gro filestream and returns topology of the system

        Args:
            f (FileIO): descriptor of opened file

        Returns:
            tuple(np.array, dict): `xyz` array of coordinates and topology of the system: atomic numbers `z`, names of the residues, 
            their location (indices across `xyz` and `z`)
        """
        n_atoms = int(f.readline().strip())
        residues = []
        res_ = ''
        indices = []
        loc = []
        xyzs = []
        zs = []
        for j in range(n_atoms):
            l = [i.strip() for i in f.readline().split(' ') if i != '']
            residue = l[0]
            res_index, res_name = re.split('(\d+)(\w+)', residue)[1:3]
            if residue != res_:
                residues.append(res_name)
                indices.append(res_index)
                res_ = residue
                loc.append(j)
            at = l[1]
            if len(at) > 3:
                xyz = [10*float(i) for i in l[2:5]]
            else:
                xyz = [10*float(i) for i in l[3:6]]
            z = _ATOMS[re.split('(\D+)(\d+)?', at)[1]]
            zs.append(z)
            xyzs.append(xyz)
        l = [i.strip() for i in f.readline().split(' ') if i != '']
        box = [10*float(i) for i in l]
        loc.append(n_atoms)
        return np.array(xyzs), {'z': np.array(zs), 'residues': residues, 'indices': indices, 'loc': np.array(loc), 'box': np.array(box)}
    
    @classmethod
    def _parse_file(cls, filename):
        with open(filename) as f:
            f.readline()
            topo = cls._parse(f)
        return topo
    
    def from_topology(self, xyz, topo):
        """ Parameters
        ----------
        xyz: numpy.ndarray
        topo: dict
            z: numpy.ndarray or list
                atomic numbers
            residues: list
                names of the residues
            indices: list
                indices of the residues
            loc: numpy.ndarray
                location of the residues within the xyz, N + 1
            box: numpy.ndarray
                box size
        box: int
            
        """
        self.box = topo['box']
        n = len(topo['indices'])
        s1 = topo['loc'][0]
        residues = topo['residues']
        indices = topo['indices']
        z = topo['z']
        for j in range(n):
            s2 = topo['loc'][j+1]
            mol = GMolecule(np.array(xyz[s1:s2]), z[s1:s2], residues[j], indices[j])
            # here potential optimization: fix only the molecules which are close to the box edge
            mol.fix(self.box)
            self.molecules.append(mol)
            s1 = s2
        self.indices = [molecule.index for molecule in self.molecules]
        self.n_molecules = len(self.molecules)
        mol_names = [molecule.name for molecule in self.molecules]
        mol_unique = set(mol_names)
        self.n_molecules_ = {name: mol_names.count(name) for name in mol_unique}
        return
    
    def from_file(self, filename):
        self.from_topology(*(self._parse_file(filename)))
    
    def from_molecules(self, molecules):
        self.molecules = molecules
        self.indices = [molecule.index for molecule in self.molecules]
        self.n_molecules = len(self.molecules)
        mol_names = [molecule.name for molecule in molecules]
        mol_unique = set(mol_names)
        self.n_molecules_ = {name: mol_names.count(name) for name in mol_unique}
        
    def __repr__(self):
        return f'GFrame object with {self.n_molecules} molecules'
    
    @property
    def xyz(self):
        return np.concatenate([molecule.xyz for molecule in self.molecules])
    
    def _xyz(self, residue):
        return np.concatenate([molecule.xyz for molecule in self.molecules if molecule.name == residue])
    
    @property
    def z(self):
        return np.concatenate([molecule.z for molecule in self.molecules])
    
    def _z(self, residue):
        return np.concatenate([molecule.z for molecule in self.molecules if molecule.name == residue])
    
    def center(self, n_molecule, n_atom):
        """Centers the structures around given molecule and atom"""
        xyz0 = self.molecules[n_molecule].xyz[n_atom, :]
        molecules = [molecule.center(xyz0, self.box) for molecule in self.molecules]
        return molecules
    
    def closest_atoms(self, idx, dist, center_at=0):
        """ For the selected molecule (`idx`), find all neighbours within a specific distance `dist`.
            The molecule is a neighbour, if at least one atom is closer than `dist` from any atom of the selected molecule.
            The entire neighbor is included into an output. Output system is centered so that idx_atom is [0.0, 0.0, 0.0]
        """
        sub = [idx]
        molecules = self.center(idx, center_at)
        xyz0 = molecules[idx].xyz
        for j in range(self.n_molecules):
            if j != idx:
                xyz1 = molecules[j].xyz
                d = cdist(xyz0, xyz1).ravel()
                if d.min() < dist:
                    sub.append(j)
        new = type(self)([molecules[i] for i in sub])
        new.box = self.box
        return new
    
    def neighbourhood(self, idx, dist: dict, center_at=0):
        """Finds all the molecules fully within a given distance `dist` from the molecule `idx`. 
        `dist` can be set different for each molecule type (via using dict).
        Output system is centered so that idx_atom is [0.0, 0.0, 0.0]
        TODO: add possibility to specify molecule by its name
        """
        sub = [idx]
        molecules = self.center(idx, center_at)
        xyz0 = molecules[idx].xyz
        for j in range(self.n_molecules):
            if j != idx:
                xyz1 = molecules[j].xyz
                molname = molecules[j].name
                d = cdist(xyz0, xyz1)
                if type(dist) == int:
                    d_threshold = dist
                else:
                    d_threshold = dist.get(molname) or 0
                if all(d.min(axis=0) < d_threshold):
                    sub.append(j)
        new = type(self)([molecules[i] for i in sub])
        new.box = self.box
        return new
    
    def to_xyz(self, filename=None, comment=''):
        """Prepare for output"""
        names = ','.join([molecule.fullname for molecule in self.molecules])
        xyz = self.xyz
        z = self.z
        n = len(z)
        lines = [str(n), '# ' + comment + ' #MOL=' + names]
        for j in range(n):
            lines.append(f'{z[j]:3d}{xyz[j, 0]:8.3f}{xyz[j, 1]:8.3f}{xyz[j, 2]:8.3f}')
        if filename is None:
            for line in lines:
                print(line)
        elif isinstance(filename, str):
            if not filename.endswith('xyz'):
                filename = filename + '.xyz'
            with open(filename, 'w') as f:
                for line in lines:
                    print(line, file=f)
    
    def to_gro(self, filename=None, comment=''):
        """Saves gro object into *gro  legacy format. Note that coordinates are in nm not angstrom. """
        n = len(self.z)
        lines  = [comment, str(n)]
        atom_idx = 0
        for molecule in self.molecules:
            molecule_idx =  molecule.index
            molecule_name = molecule.name
            xyz = molecule.xyz
            z = molecule.z
            for j in range(len(z)):
                atom_idx += 1
                atom_name = f'{_CHARGES[z[j]]}{j+1}'
                xyz_ = xyz[j]
                s = f'{molecule_idx:>5}{molecule_name:<5}{atom_name:>5}{atom_idx:>5}{0.1*xyz_[0]:8.3f}{0.1*xyz_[1]:8.3f}{0.1*xyz_[2]:8.3f}'
                lines.append(s)
        box = 10 if self.box is None else self.box[0] 
        lines.append(f'{0.1*box:8.3f}{0.1*box:8.3f}{0.1*box:8.3f}')
        if filename is None:
            for line in lines:
                print(line)
            return
        with open(filename, 'w') as f:
            for line in lines:
                print(line, file=f)
        
    
    def fix(self):
        """TODO: call fix() for each molecule"""
        pass
    
    
class GTrajectory:
    def __init__(self, **kwargs):
        self.snapshots = []
        self.t = []
        if kwargs.get('file_gro'):
            if kwargs.get('file_trr'):
                self.from_trr(file_trr=kwargs['file_trr'], file_gro=kwargs['file_gro'])
            else:
                self.from_gro(file_gro=kwargs['file_gro'])

    def from_gro(self, file_gro):
        self.snapshots = []
        self.t = []
        with open(file_gro) as f:
            l = f.readline()
            while l != '':
                if 't=' in l:
                    s = l.index('t=')
                    p = [i for i in l[s+2:].split(' ') if i != '']
                    t = float(p[0])
                else:
                    t = None
                gro = GFrame()
                gro.t = t
                xyz, topo = GFrame._parse(f)
                gro.from_topology(xyz, topo)
                self.snapshots.append(gro)
                l = f.readline()
                self.t.append(t)
    
    def __getitem__(self, idx):
        return self.snapshots[idx]
    
    def __repr__(self):
        return f'Trajectory with {self.n_frames} snapshots'
    
    @property
    def n_frames(self):
        return len(self.snapshots)
    
    def from_trr(self, file_trr, file_gro):
        # trr file is just binary data file, 120 bytes os header + 4 x 3 x (N_atoms) bytes, big-endian (>f)
        # read gro, get atoms -> load corresponding xyz
        self.snapshots = []
        self.t = []
        xyz, topo = GFrame._parse_file(file_gro)
        with open(file_trr, 'rb') as f:
            h = f.read(120)
            while len(h) == 120:
                N,  = struct.unpack('>I', h[52:56])
                Na, frame = struct.unpack('>2I', h[64:72])
                t,  = struct.unpack('>f', h[76:80])
                box = struct.unpack('>9f', h[84:120])[::4]
                topo['box'] = 10*np.array(box)  # need to be set in advance for correct fix of the molecules
                assert N == 12*Na, 'Some mismatch in binary file'
                d = f.read(N)
                xyz = 10*np.array(struct.unpack(f'>{3*Na}f', d)).reshape(Na, 3)
                gro = GFrame()
                gro.from_topology(xyz, topo)
                gro.t = t
                self.snapshots.append(gro)
                self.t.append(t)
                h = f.read(120)
