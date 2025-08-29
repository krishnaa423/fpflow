#region modules
from fpflow.structure.struct import Struct
import numpy as np 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class WannierStruct:
    def __init__(self, struct: Struct, struct_idx: int = 0):
        self.struct: Struct = struct
        self.struct_idx: int = struct_idx

    @property
    def cell(self):
        return {'unit': 'Ang', 'data': self.struct.atoms[self.struct_idx].get_cell().array.tolist()}

    @property
    def atomic_positions(self):
        ang_coords = self.struct.atoms[self.struct_idx].get_positions()
        cell = self.struct.atoms[self.struct_idx].get_cell().array
        data = (ang_coords @ np.linalg.inv(cell)).tolist()
        return {'data': data}

#endregion