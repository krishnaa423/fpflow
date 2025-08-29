#region modules
from fpflow.structure.struct import Struct
from ase.data import chemical_symbols, atomic_masses
import numpy as np 
import os
from importlib.util import find_spec 
from lxml import etree 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class AbacusStruct(Struct):
    @property
    def atomic_species(self):
        data = []
        
        for atm_num in np.unique(self.atoms[self.struct_idx].get_atomic_numbers()):
            data.append([
                chemical_symbols[atm_num],
                float(atomic_masses[atm_num]),
                f'{chemical_symbols[atm_num]}.upf'
            ])

        return data 

    @property
    def orbitals(self):
        output = []
        
        for atm_num in np.unique(self.atoms[self.struct_idx].get_atomic_numbers()):
            output.append([
                f'{chemical_symbols[atm_num]}.orb'
            ])

        return output 

    @property
    def cell(self):
        return self.atoms[self.struct_idx].get_cell().array.tolist()

    @property
    def atomic_positions(self):
        data = []

        pos_array = self.atoms[self.struct_idx].get_positions()
        numbers = self.atoms[self.struct_idx].get_atomic_numbers()

        for atm_num, row in zip(numbers, pos_array):
                data.append([
                    chemical_symbols[atm_num],
                    float(row[0]),
                    float(row[1]),
                    float(row[2]),
                    0.0, # move_x,
                    0.0, # move_y,
                    0.0, # move_z,
                    0.0, # magnetic moment.
                ])
        
        return {'unit': 'Cartesian', 'data': data}
    

    @property
    def lattice_constant(self):
        return [[1.0]]
    
    def get_pseudos_list(self, xc: str=None, is_soc: bool=False):
        string_list = []

        pkg_dir = os.path.dirname(find_spec('fpflow').origin)
        string_list.append(pkg_dir)
        string_list.append('data')
        string_list.append('pseudos')
        string_list.append('abacus')
        string_list.append('')  # Placeholder for the element symbol. 

        paths = []
        filenames = []
        for symbol in self.atoms[self.struct_idx].symbols:
            string_list[-1] = f'{symbol}.upf'
            path = os.path.join(*string_list)
            filename = os.path.basename(path)
            paths.append(path)
            filenames.append(filename)

        return paths, filenames
    
    def max_val(self, xc: str=None, is_soc: bool=False):
        paths, filenames = self.get_pseudos_list(xc=xc, is_soc=is_soc)
    
        total_electrons: int = 0
        for path in paths:
            # Add number of electrons from each atom type. 
            total_electrons += int(float(etree.parse(path).xpath('string(//@z_valence)').strip()))

        return total_electrons

#endregion