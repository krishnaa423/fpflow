#region modules
from ase import Atoms 
from ase.io import read, write
from ase.units import Angstrom, Bohr
from ase.data import atomic_numbers, chemical_symbols
from fpflow.inputs.inputyaml import InputYaml
import jmespath
from typing import List 
import numpy as np 
import operator 
from fpflow.io.logging import get_logger
from importlib.util import find_spec
import os 
from fpflow.structure.external_struct.mp import MateralsProject
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class CellArray:
    def __init__(self, array: np.array):
        self.array: np.ndarray = array 

    @classmethod
    def from_structdict(cls, dict: dict):
        assert operator.xor(jmespath.search('cell.vectors', dict)  is not None, jmespath.search('cell.bravais_lattice_info', dict)  is not None), 'vectors or bravais lattice only should be present.'

        read_if_relaxed: bool = False if jmespath.search('cell.read_if_relaxed', dict) is None else jmespath.search('cell.read_if_relaxed', dict)

        # Read relaxed file if present. 
        if read_if_relaxed and os.path.exists('./relaxed_cell_parameters.txt'):
            vectors = np.loadtxt('./relaxed_cell_parameters.txt')
            return cls(vectors)

        unit: str = jmespath.search('cell.unit', dict)
        vectors: np.array = None if jmespath.search('cell.vectors[*]', dict) is None else np.array(jmespath.search('cell.vectors[*]', dict))
        ibrav = jmespath.search('cell.bravais_lattice_info.ibrav', dict)
        A = jmespath.search('cell.bravais_lattice_info.A', dict)
        B = jmespath.search('cell.bravais_lattice_info.B', dict)
        C = jmespath.search('cell.bravais_lattice_info.C', dict)

        factor: float = None 
        match unit:
            case 'angstrom': 
                factor = Angstrom
            case 'bohr':
                factor = Bohr
            case 'alat':
                factor = A

        if vectors is not None:
            return cls(vectors*factor)
        
        # If no cell vectors, use bravais lattice info. 
        match ibrav:
            case 'sc':
                vectors = np.diag([A, A, A])
            case 'fcc':
                vectors = np.array([
                    [-0.5, 0.0, 0.5],
                    [0.0, 0.5, 0.5],
                    [-0.5, 0.5, 0.0],
                ])*A
            case 'bcc':
                vectors = np.array([
                    [0.5, 0.5, 0.5],
                    [-0.5, 0.5, 0.5],
                    [-0.5, -0.5, 0.5],
                ])*A
            case 'tetra':
                vectors = np.diag([A, A, C])
            case 'ortho':
                vectors = np.diag([A, B, C])
            case 'hex':
                vectors = np.array([
                    [A, 0.0, 0.0],
                    [-0.5*A, 0.5*np.sqrt(3)*A, 0.0],
                    [0.0, 0.0, C],
                ])
            case _:
                vectors = np.diag([1.0, 1.0, 1.0])
        return CellArray(vectors)

class PositionArray:
    def __init__(self, array: np.array):
        self.array: np.ndarray = array 

    @classmethod
    def from_structdict(cls, dict: dict):
        cell_array = CellArray.from_structdict(dict)
        
        # Read relaxed file if present. 
        read_if_relaxed: bool = False if jmespath.search('positions.read_if_relaxed', dict) is None else jmespath.search('positions.read_if_relaxed', dict)
        if read_if_relaxed and os.path.exists('./relaxed_atomic_positions.txt'):
            vectors = np.loadtxt('./relaxed_atomic_positions.txt', usecols=(1, 2, 3))
            return cls(vectors)

        unit: str = jmespath.search('positions.unit', dict)
        coords: np.ndarray = np.array(jmespath.search('positions.coords[*][1:4]', dict))
        A: float = jmespath.search('cell.bravais_lattice_info.A', dict)

        match unit:
            case 'angstrom': 
                coords = coords * Angstrom
            case 'bohr':
                coords = coords * Bohr
            case 'alat':
                coords = coords * A
            case 'crystal':
                coords = coords @ cell_array.array

        return cls(coords)

class AtomicNumberArray:
    def __init__(self, array: np.array):
        self.array: np.ndarray = array 

    @classmethod
    def from_structdict(cls, dict: dict):
        symbols = jmespath.search('positions.coords[*][0]', dict)
        numbers: np.ndarray = np.array([atomic_numbers[s] for s in symbols])
        return AtomicNumberArray(numbers)

    def get_atomic_symbols(self) -> List[str]:
        symbols = [chemical_symbols[Z] for Z in self.array]
        return symbols

class Struct:
    def __init__(self, atoms: List[Atoms], struct_idx:int=0):
        self.atoms: List[Atoms] = atoms 
        self.struct_idx: int = struct_idx

    @classmethod
    def from_yaml_file(cls, filename: str, **kwargs):
        inputdict: dict = InputYaml.from_yaml_file(filename).inputdict
        atoms: List[Atoms] = cls.get_atoms_from_inputdict(inputdict)
        return cls(atoms=atoms, **kwargs)
    
    @classmethod
    def from_inputdict(cls, inputdict: dict, **kwargs):
        atoms: List[Atoms] = cls.get_atoms_from_inputdict(inputdict)
        struct_idx: int = jmespath.search('structures.active_idx', inputdict)
        return cls(atoms=atoms, struct_idx=struct_idx, **kwargs)
    
    @classmethod
    def _get_external_atoms(cls, struct_dict: dict) -> Atoms:
        struct_type: str = jmespath.search('file.type', struct_dict)
        struct_value: str = jmespath.search('file.value', struct_dict)

        match struct_type:
            case 'localfile':
                return read(struct_value)
            case 'mpapi-id':
                return MateralsProject.from_id(struct_value).atoms
            case 'mpapi-formula':
                return MateralsProject.from_formula(struct_value).atoms

    @classmethod
    def get_atoms_from_inputdict(cls, inputdict: dict) -> List[Atoms]:
        structures = jmespath.search('structures.list[*]', inputdict)
        atoms: List[Atoms] = []

        for struct_idx, struct_item in enumerate(structures):
            # If file is present. 
            if struct_item['file']  is not None:
                structure: Atoms = cls._get_external_atoms(struct_dict=struct_item)
            else: #otherwise. 
                structure: Atoms = Atoms(
                    numbers=AtomicNumberArray.from_structdict(struct_item).array,
                    positions=PositionArray.from_structdict(struct_item).array,
                    cell=CellArray.from_structdict(struct_item).array,
                    pbc=[True, True, True],
                )
            
            # Create supercell if specified. 
            supercell_size: list = jmespath.search(f'structures.list[{struct_idx}].supercell_size', inputdict)
            if supercell_size is not None:
                structure = structure * supercell_size
            
            # Append the structure. 
            atoms.append(structure)

            # Write structure file. 
            write(f'atoms_{struct_idx}.xsf', structure)

        return atoms
    
    
    def ntyp(self, idx:int=0):
        return len(np.unique(self.atoms[idx].get_atomic_numbers()))

    
    def nat(self, idx:int=0):
        return len(self.atoms[idx]) 

#endregion