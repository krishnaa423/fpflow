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
    def from_special_pts_list(cls, atoms: Atoms, special_points: Iterable[str], npoint_per_segment: int = 20):
        special_points_loc = get_special_points(atoms.cell)

        num_special_points = len(special_points)
        kpts = np.zeros(shape=((num_special_points-1)*npoint_per_segment+1, 3), dtype='f8')

        # Add points between the special points. 
        for sp_idx in range(num_special_points-1):
            for coord in range(3):
                start = special_points_loc[special_points[sp_idx]][coord]
                stop = special_points_loc[special_points[sp_idx+1]][coord]
                step = (stop - start)/npoint_per_segment
                kpts[sp_idx*npoint_per_segment:(sp_idx+1)*npoint_per_segment, coord] = np.arange(start, stop, step) if step!=0.0 else 0.0

        # Add the final kpoint. 
        kpts[-1, :] = np.array(special_points_loc[special_points[-1]])

        return cls(
            special_points=special_points,
            npoints_segment=npoint_per_segment,
            kpts=kpts,
            atoms=atoms,
        )

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
    def ase_axis(self):
        bandpath: BandPath = self.atoms.cell.bandpath(
            path=''.join(self.special_points), 
            npoints=(len(self.special_points)-1)*self.npoints_segment+1
        )

        return bandpath.get_linear_kpoint_axis()
    
    @property
    def even_spaced_axis(self):
        nseg = len(self.special_points) - 1
        n = nseg * self.npoints_segment + 1
        xaxis = np.arange(n, dtype=float)
        
        # xticks at the indices of each special point
        xticks = [i * self.npoints_segment for i in range(nseg)] + [n - 1]
        xlabels = list(self.special_points)
        
        return xaxis, xticks, xlabels
    
    @property
    def phonopy_axis(self):
        nseg = len(self.special_points) - 1
        n = nseg * self.npoints_segment  # e.g., 2 * 10 = 20

        # Build cumulative path without duplicating junctions:
        # - First segment: include left endpoint (L), exclude right (Γ)
        # - Subsequent segments: exclude left, include right (last includes X)
        parts = []
        for i in range(nseg):
            left, right = i, i + 1
            if i == 0:
                seg = np.linspace(left, right, self.npoints_segment, endpoint=False)            # 0.0..0.9
            else:
                seg = np.linspace(left, right, self.npoints_segment + 1, endpoint=True)[1:]     # 1.1..2.0
            parts.append(seg)

        xaxis = np.concatenate(parts).astype(float)

        xticks = np.arange(nseg + 1, dtype=float)         # [0., 1., 2.]
        xlabels = list(self.special_points)               # ["L", "G", "X"]
        return xaxis, xticks, xlabels

    @property
    def matdyn_str(self):
        output = ''
        special_points = get_special_points(self.atoms.cell)

        output += f'{len(self.special_points)}\n'

        for path_special_point in self.special_points:
            coord = special_points[path_special_point]
            output += f'{coord[0]:15.10f} {coord[1]:15.10f} {coord[2]:15.10f} {self.npoints_segment} !{path_special_point}\n'
        
        return output 

    @property
    def dftelbands_list(self):
        output = []
        special_points = get_special_points(self.atoms.cell)

        for path_special_point in self.special_points:
            coord = special_points[path_special_point]
            output.append([
                coord[0],
                coord[1],
                coord[2],
                self.npoints_segment,
            ])

        return output
        
    @property
    def wannierpath_list(self):
        output = []

        for sp_idx in range(len(self.special_points)-1):
            from_sp = self.special_points[sp_idx]
            to_sp = self.special_points[sp_idx+1]
            from_coord = get_special_points(self.atoms.cell)[from_sp]
            to_coord = get_special_points(self.atoms.cell)[to_sp]
            output.append([
                from_sp,
                from_coord[0],
                from_coord[1],
                from_coord[2],
                to_sp,
                to_coord[0],
                to_coord[1],
                to_coord[2],
            ])

        return output


#endregion