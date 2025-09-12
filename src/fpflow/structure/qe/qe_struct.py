#region modules
from fpflow.structure.struct import Struct
import numpy as np
from ase.data import chemical_symbols, atomic_masses
import os
from importlib.util import find_spec 
from lxml import etree 
import jmespath
from fpflow.io.logging import get_logger
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class QeStruct(Struct):
    @property
    def cell(self):
        return {'unit': 'angstrom', 'data': self.atoms[self.struct_idx].get_cell().array.tolist()}

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
    def atomic_positions(self):
        numbers = self.atoms[self.struct_idx].get_atomic_numbers()
        data = self.atoms[self.struct_idx].get_positions().tolist()

        for symbol, pos in zip(numbers, data):
            pos.insert(0, chemical_symbols[symbol])

        return {'unit': 'angstrom', 'data': data}
    
    def get_pseudos_list(self, xc: str=None, is_soc: bool=False):
        string_list = []

        pkg_dir = os.path.dirname(find_spec('fpflow').origin)
        string_list.append(pkg_dir)
        string_list.append('data')
        string_list.append('pseudos')
        string_list.append('qe')
        
        sub_string = ''
        sub_string += 'fr' if is_soc else 'sr'
        sub_string += '_'
        sub_string += f'{xc}' if xc is not None else 'pbe'

        string_list.append(sub_string)

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
    
        total_val_bands: int = 0
        for path in paths:
            # Add number of electrons from each atom type. 
            total_val_bands += int(float(etree.parse(path).xpath('string(//@z_valence)').strip()))

        if not is_soc:
            total_val_bands /= 2

        return total_val_bands
    
    def get_ibrav(self, inputdict: dict) -> int:
        struct_dict: dict = jmespath.search('structures.list[*]', inputdict)[self.struct_idx]
        ibrav_str: str = 'free' if jmespath.search('cell.bravais_lattice_info.ibrav', struct_dict) is None else jmespath.search('cell.bravais_lattice_info.ibrav', struct_dict)

        match ibrav_str:
            case 'sc':
                return 1
            case 'fcc':
                return 2
            case 'bcc':
                return 3
            case 'tetra':
                return 4
            case 'ortho':
                return 8
            case 'hex':
                return 4
            case _:
                return 0

#endregion