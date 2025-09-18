#region modules.
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.nestedstep_base import NestedBaseStep
import jmespath 
from fpflow.schedulers.scheduler import Scheduler
import copy 
import numpy as np 
from glom import Path 
from fpflow.plots.phbands import PhonopyPlot
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms 
from ase import Atoms
from fpflow.structure.struct import Struct

#endregion

#region variables.
#endregion

#region functions.
#endregion

#region classes.
class QePhonopyNestedStep(NestedBaseStep):
    def __init__(self, **kwargs):
        super().__init__(steptag='phonopy', **kwargs)

    def get_phonopy_dispatoms(self):
        # Get phonopy atoms. 
        struct: Struct = Struct.from_inputdict(self.inputdict)
        ase_atoms: Atoms = struct.atoms[struct.struct_idx]
        phonopy_atoms: PhonopyAtoms = PhonopyAtoms(
            symbols=ase_atoms.get_chemical_symbols(),
            cell=ase_atoms.get_cell(),
            scaled_positions=ase_atoms.get_scaled_positions()
        )
        supercell_size = jmespath.search('dfpt.qgrid', self.inputdict)

        # Generate disps. 
        phonon = Phonopy(phonopy_atoms, supercell_matrix=supercell_size)
        phonon.generate_displacements(distance=0.01)
        phonon.save('phonopy_disp.yaml')

        # output dict. 
        dispatoms_dict: dict = {}
        for idx, sc in enumerate(phonon.supercells_with_displacements):
            dispatoms_dict[idx] = sc

        return dispatoms_dict
    
    @property
    def folder_inputdict_changes_map(self) -> dict:
        '''
        Return:
            {
                'folder1': [
                    {
                        'path': 'some.path1',
                        'value': 12,
                    },
                    {
                        'path': 'some.path2',
                        'value': 45,
                    },
                ],
                'folder2': [
                    {
                        'path': 'some.path1',
                        'value': 123,
                    },
                    {
                        'path': 'some.path2',
                        'value': 345,
                    },
                ],
            }
        '''
        output: dict = {}
        dispstructs_dict: dict = self.get_phonopy_dispatoms()

        for key, sc in dispstructs_dict.items():
            symbol_and_positions = []
            for symbol, pos in zip(sc.get_chemical_symbols(), sc.get_scaled_positions()):
                symbol_and_positions.append([
                    symbol,
                    float(pos[0]),
                    float(pos[1]),
                    float(pos[2]),
                ])

            output[f'./phonopy/dset_{key}'] = [
                {
                    'path':  'generator.pre_steps',
                    'value': ['pseudos_qe'],
                },
                {
                    'path':  'generator.steps',
                    'value': ['scf_qe'],
                },
                {
                    'path': 'manager.steps',
                    'value': ['scf_qe'],
                },
                {
                    'path': 'manager.plots',
                    'value': ['scf_qe'],
                },
                {
                    'path': 'structures.active_idx',
                    'value': 0,
                },
                {
                    'path': 'structures.list[0].file',
                    'value': None,
                },
                {
                    'path': 'structures.list[0].cell.vectors',
                    'value': sc.cell.tolist(),
                },
                {
                    'path': 'structures.list[0].cell.bravais_lattice_info',
                    'value': None,
                },
                {
                    'path': 'structures.list[0].positions.coords',
                    'value': symbol_and_positions,
                },
                {
                    'path': 'structures.list[0].supercell_size',
                    'value': None,
                },
                {
                    'path': 'scf.kgrid',
                    'value': [1, 1, 1],
                },
                {
                    'path': 'scf.args.system.nosym',
                    'value': True,
                },
                {
                    'path': 'scf.job_info',
                    'value': jmespath.search('phonopy.job_info', self.inputdict),
                },
            ]

            return output 
    
    @property
    def extra_filecontents(self) -> dict:
        output: dict = super().extra_filecontents

        return output
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_phonopy.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './esr',
            './job_phonopy.sh',
            './script_phonopy.py',
            './script_phonopy.py.out',
        ]

#endregion
