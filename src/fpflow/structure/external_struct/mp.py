#region modules
import os 
from ase import Atoms 
from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
#endregion

#region variables
#endregion

#region functions
mpid = os.environ['MP_API']
#endregion

#region classes
class MateralsProject:
    def __init__(self, **kwargs):
        self.atoms: Atoms = None 

        for key, value in kwargs.items():
            setattr(self, key, value)


    @classmethod
    def from_id(cls, id: int):
        with MPRester(mpid) as mpr:
            docs = mpr.materials.summary.search(material_ids=[id])
            atoms: Atoms = AseAtomsAdaptor().get_atoms(docs[0].structure)

        return cls(atoms=atoms)
    
    @classmethod
    def from_formula(cls, formula: int):
        with MPRester(mpid) as mpr:
            docs = mpr.materials.summary.search(formula=formula, energy_above_hull=[0, 0])
            atoms: Atoms = AseAtomsAdaptor().get_atoms(docs[0].structure)

        return cls(atoms=atoms)

#endregion