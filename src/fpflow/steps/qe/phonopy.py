#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from importlib.util import find_spec
import glob 
from ase.dft.kpoints import get_special_points
from fpflow.structure.struct import Struct
from ase import Atoms
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.plots.phbands import PhonopyPlot
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms 
import copy 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePhonopyStep(Step):
    @property
    def script_phonopy_fc_bs(self):
        regex_pattern: str = r"phonopy-supercell-(\d+)\.out$"

        return f'''
# script_phonopy_fc_bs.py (patched)
from fpflow.inputs.inputyaml import InputYaml
from phonopy import load
from phonopy.interface.qe import parse_set_of_forces
from ase.dft.kpoints import get_special_points
import numpy as np
import jmespath
import glob
import re
from fpflow.structure.struct import Struct
from ase import Atoms
import os

def sorted_supercell_outs(pattern="phonopy-supercell-*.out"):
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No QE output files matched pattern")
    # sort by trailing integer index: phonopy-supercell-<idx>.out
    def idx_of(f):
        m = re.search(r"{regex_pattern}", f)
        return int(m.group(1)) if m else 1_000_000
    return sorted(files, key=idx_of)

# 1) Load displacements + structure that Phonopy created earlier
phonon = load("phonopy_disp.yaml")  # has supercell, displacements, etc.

# 2) Collect ONLY QE .out files in the correct order and get num_atoms
qe_outs = sorted_supercell_outs("phonopy-supercell-*.out")
num_atoms = len(phonon.supercell)

# 3) Parse forces from QE outputs
#    Your Phonopy version expects (forces_filenames, num_atoms, ...)
forces = parse_set_of_forces(
    forces_filenames=qe_outs,
    num_atoms=num_atoms
)

# 4) Attach forces and build force constants
phonon.forces = forces
phonon.produce_force_constants()

# 5) Save for later reuse
phonon.save("force_constants.h5")
phonon.save("phonopy.yaml")

# 6) Build band path using ASE special k-points and the kpath from input yaml
inputdict: dict = InputYaml.from_yaml_file(yaml_filename='./input.yaml').inputdict
struct: Struct = Struct.from_inputdict(inputdict=inputdict)
ase_atoms: Atoms = struct.atoms[struct.struct_idx]
points = get_special_points(ase_atoms.cell)
path = jmespath.search('kpath.special_points', inputdict)  # e.g., ["G","X","W",...]
nseg = jmespath.search('kpath.npoints_segment', inputdict)  # integer

bands = [[points[s] + (points[e] - points[s]) * i / nseg
          for i in range(nseg + 1)]
         for s, e in zip(path[:-1], path[1:])]

# 7) Compute & write band structure
phonon.run_band_structure(bands)
phonon.write_yaml_band_structure("band.yaml")
'''
        
    @property
    def qedisp_files(self) -> list[str]:
        qedisp_files_dict: dict = {}

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
        
        # Get supercells.
        for idx, sc in enumerate(phonon.supercells_with_displacements):
            ase_atoms = Atoms(
                numbers=sc.get_atomic_numbers(),
                positions=sc.get_positions(),
                cell=sc.cell,
                pbc=[True, True, True],
            )
            qestruct: QeStruct = QeStruct(atoms=[ase_atoms])

            qedict: dict = {
                'control': {
                    'outdir': './tmp',
                    'prefix': 'struct',
                    'pseudo_dir': './pseudos/qe',
                    'calculation': 'scf',
                    'tprnfor': True,
                },
                'system': {
                    'ibrav': 0,
                    'ntyp': qestruct.ntyp(),
                    'nat': qestruct.nat(),
                    'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                    'nosym': True,
                },
                'electrons': {},
                'ions': {},
                'cell': {},
                'atomic_species': qestruct.atomic_species,
                'cell_parameters': qestruct.cell,
                'atomic_positions': qestruct.atomic_positions,
                'k_points': {
                    'type': 'automatic',
                    'data': [
                        1,
                        1,
                        1,
                        0,
                        0,
                        0,
                    ],
                }
            }
            if jmespath.search('scf.is_spinorbit', self.inputdict):
                qedict['system']['noncolin'] = True
                qedict['system']['lspinorb'] = True

            # Update if needed. 
            update_dict(qedict, jmespath.search('scf.args', self.inputdict))

            qedisp_files_dict[f'phonopy-supercell-{idx}'] = QeGrammar().write(qedict)

        return qedisp_files_dict 

    @property
    def job_phonopy(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dfpt.job_info')

        qefiles: list[str] = self.qedisp_files.keys()

        # Create supercell job.
        start_idx = 0
        stop_idx = len(qefiles)
        debug_str: str = '\n'
        files_bashvar_str: str = '\nfiles=('
        files_args_str: str = ''
        for file_idx, file in enumerate(qefiles): 
            files_bashvar_str += f'"{file}" '
            files_args_str += f' {file}.out '
            debug_str += f'#idx: {file_idx}, filename: {file}\n'
        files_bashvar_str += ')\n\n'
        debug_str += '\n\n'
        file_variable = '${files[$i]}'

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{debug_str}

{files_bashvar_str}

start={start_idx}
stop={stop_idx-1}

for (( i=$start; i<=$stop; i++ )); do
{scheduler.get_exec_prefix()}pw.x < {file_variable} &> {file_variable}.out
done

# Post processing. This should create force_constants.yaml and band.yaml
python script_phonopy_fc_bs.py &> script_phonopy_fc_bs.py.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        output: dict = {
            'job_phonopy.sh': self.job_phonopy,
            'script_phonopy_fc_bs.py': self.script_phonopy_fc_bs,
        }
        output.update(self.qedisp_files)

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
            './phonopy_disp.yaml',
            './phonopy-supercell*',
            './job_phonopy.sh',
            './phonopy.yaml',
            './force_constants.h5',
            './band.yaml',
            './script_phonopy_fc_bs.py',
            './script_phonopy_fc_bs.py.out',
        ]
    
    def plot(self):
        PhonopyPlot().save_figures()

#endregion