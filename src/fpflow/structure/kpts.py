#region modules
import numpy as np 
from typing import Iterable
import itertools 
from fpflow.structure.struct import Struct
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.io.logging import get_logger
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
            case True:
                class_call = Kpts
            case False:
                class_call = FbzKpts

        result = class_call(
            kgrid=kgrid,
            qshift=qshift,
            is_reduced=is_reduced,
            is_grid=True,
        ) 
        result.populate_kpts(**kwargs)

        return result 

class FbzKpts(Kpts):
    def populate_kpts(self, **kwargs):
        raise NotImplementedError()
    
#endregion