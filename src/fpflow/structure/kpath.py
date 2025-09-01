#region modules
from typing import List, Iterable
from fpflow.structure.struct import Struct
from fpflow.inputs.inputyaml import InputYaml
import jmespath
from ase import Atoms 
from ase.dft.kpoints import BandPath, get_special_points
import numpy as np 
from fpflow.io.logging import get_logger
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class Kpath:
    def __init__(self, **kwargs):
        self.special_points: Iterable[str] = None 
        self.npoints_segment: int = None 
        self.kpts: Iterable = None 
        self.atoms: Atoms = None 

        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def from_special_pts_list(cls, special_points: Iterable[str], npoint_per_segment: int = 20):
        raise NotImplementedError()

    @classmethod
    def from_yamlfile(cls, filename: str='./input.yaml', struct_idx: int =0):
        struct: Struct = Struct.from_yaml_file(filename)
        struct_idx: int = struct_idx
        atoms: Atoms = struct.atoms[struct_idx]
        inputdict: dict = InputYaml.from_yaml_file(filename).inputdict

        special_points: Iterable[str] = jmespath.search('kpath.special_points[*]', inputdict)
        npoints_segment: int = jmespath.search('kpath.npoints_segment', inputdict)

        special_points_loc = get_special_points(atoms.cell)

        num_special_points = len(special_points)
        kpts = np.zeros(shape=((num_special_points-1)*npoints_segment+1, 3), dtype='f8')

        # Add points between the special points. 
        for sp_idx in range(num_special_points-1):
            for coord in range(3):
                start = special_points_loc[special_points[sp_idx]][coord]
                stop = special_points_loc[special_points[sp_idx+1]][coord]
                step = (stop - start)/npoints_segment
                kpts[sp_idx*npoints_segment:(sp_idx+1)*npoints_segment, coord] = np.arange(start, stop, step) if step!=0.0 else 0.0

        # Add the final kpoint. 
        kpts[-1, :] = np.array(special_points_loc[special_points[-1]])

        return cls(
            special_points=special_points,
            npoints_segment=npoints_segment,
            kpts=kpts,
            atoms=atoms,
        )

    @property
    def axis(self):
        bandpath: BandPath = self.atoms.cell.bandpath(
            path=''.join(self.special_points), 
            npoints=len(self.special_points)*self.npoints_segment
        )

        return bandpath.get_linear_kpoint_axis()
    
    @property
    def matdyn_str(self):
        output = ''
        special_points = get_special_points(self.atoms.cell)

        output += f'{len(self.special_points)}\n'

        for path_special_point in self.special_points:
            coord = special_points[path_special_point]
            output += f'{coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f} {self.npoints_segment} !{path_special_point}\n'
        
        return output 

        
#endregion