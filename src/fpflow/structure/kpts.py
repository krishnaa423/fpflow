#region modules
import numpy as np 
from typing import Iterable
import itertools 
from fpflow.structure.struct import Struct
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.io.logging import get_logger
from ase import Atoms 
from fpflow.structure.struct import Struct
from fpflow.structure.symmetry.kgridx import get_ibz_kpts
import copy
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class Kpts:
    def __init__(self, *args, **kwargs):
        self.is_grid: bool = False
        self.is_reduced: bool = False
        self.kpts: Iterable = None 
        self.kgrid: Iterable = None 
        self.qshift: Iterable = None 

        for key, value in kwargs.items():
            setattr(self, key, value)

    def populate_kpts(self, **kwargs):
        self.kpts = np.array(list(itertools.product(range(self.kgrid[0]), range(self.kgrid[1]), range(self.kgrid[2])))) / np.array(self.kgrid)

    @classmethod
    def from_kgrid(cls, kgrid: Iterable, qshift: Iterable=[0.0, 0.0, 0.0], is_reduced: bool = False, **kwargs):
        class_call = None 
        
        match is_reduced:
            case False:
                class_call = Kpts
            case True:
                class_call = IbzKpts

        result = class_call(
            kgrid=kgrid,
            qshift=qshift,
            is_reduced=is_reduced,
            is_grid=True,
        ) 
        result.populate_kpts(**kwargs)

        return result 
    
    @property
    def nkpt(self):
        if isinstance(self.kpts, list):
            return len(self.kpts)
        elif isinstance(self.kpts, np.ndarray):
            return self.kpts.shape[0]
        else:
            return 0

    @property
    def wfn_kpts(self):
        kpts = self.kpts.tolist()
        if len(kpts[0])==3:
            for row in kpts: row.append(1.0)

        return kpts
    
    @property
    def wfnq_kpts(self):
        kpts = self.kpts.tolist()
        if len(kpts[0])==3:
            for row in kpts: row.append(1.0)

        if not isinstance(self, IbzKpts):
            for row in kpts:
                row[0] += self.qshift[0]
                row[1] += self.qshift[1]
                row[2] += self.qshift[2]

        return kpts

    @property
    def epsilon_kpts(self):
        qshift_copy = copy.deepcopy(self.qshift)
        self.qshift = [0.0, 0.0, 0.0]
        self.populate_kpts()

        kpts = self.kpts.tolist()
        if len(kpts[0])==3:
            for row in kpts: row.append(1.0); row.append(0)
        elif len(kpts[0])==4:
            for row in kpts: row[3] = 1.0; row.append(0)

        kpts[0][0] = qshift_copy[0]
        kpts[0][1] = qshift_copy[1]
        kpts[0][2] = qshift_copy[2]
        kpts[0][4] = 1

        return kpts

    @property
    def sigma_kpts(self):
        kpts = self.kpts.tolist()
        if len(kpts[0])==3:
            for row in kpts: row.append(1.0)
        elif len(kpts[0])==4:
            for row in kpts: row[3] = 1.0

        return kpts 
    
    @property
    def bseq_qpts(self):
        kpts = self.kpts.tolist()

        return kpts

class IbzKpts(Kpts):
    def populate_kpts(self, **kwargs):
        struct: Struct = Struct.from_yaml_file('./input.yaml')
        atoms: Atoms = struct.atoms[struct.struct_idx]

        # Get symmetry. 
        #TODO: Try to use spglib itself and populate the symmetry module. 
        # For now using kgrid.x 
        self.kpts = get_ibz_kpts(
            atoms=atoms,
            kgrid=self.kgrid,
            qshift=self.qshift,
        )

#endregion