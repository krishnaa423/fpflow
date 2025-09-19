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

    @property
    def script_steptag(self):
        return f'''import os
import os
import os
import sys
import glob
import subprocess
from fpflow.managers.run import subprocess_run

# # Set the stdout and stderr. 
# outfile = open('script_phonopy.py.out', 'w')
# sys.stdout = outfile
# sys.stderr = outfile

# Get the directories. 
phonopy_dirs = [inode for inode in glob.glob('./phonopy/*') if os.path.isdir(inode)]
phonopy_dirs.sort()
start: int = 0
stop: int = len(phonopy_dirs)

# Override if needed. Comment this out and set. 
#start = 0 
#stop = 1

total_time: float = 0.0
for phonopy_dir in phonopy_dirs[start:stop]:
    total_time = subprocess_run('./run.sh', total_time=total_time, dest_dir=phonopy_dir)

print(f'Done phonopy in total time: ', total_time, ' seconds.', flush=True)

# No phonopy bands stuff. 
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
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

# 1) Load displacements + structure that Phonopy created earlier
phonon = load("./phonopy_disp.yaml")  # has supercell, displacements, etc.

# 2) Collect ONLY QE .out files in the correct order and get num_atoms
qe_outs = glob.glob("./phonopy/dset_*/scf.in.out")
qe_outs.sort()
num_atoms = len(phonon.supercell)

# 3) Parse forces from QE outputs
#    Your Phonopy version expects (forces_filenames, num_atoms, ...)
forces = parse_set_of_forces(
    forces_filenames=qe_outs,
    num_atoms=num_atoms
)
print(forces)

# 4) Attach forces and build force constants
phonon.forces = forces
phonon.produce_force_constants()

# 5) Save for later reuse
phonon.save("./force_constants.h5")
phonon.save("./phonopy.yaml")

# 6) Build band path using ASE special k-points and the kpath from input yaml
inputdict: dict = InputYaml.from_yaml_file(yaml_filename='./input.yaml').inputdict
struct: Struct = Struct.from_inputdict(inputdict=inputdict)
ase_atoms: Atoms = struct.atoms[struct.struct_idx]
special_points = get_special_points(ase_atoms.cell)

# Labels & segments from your input.yaml
labels_flat = jmespath.search('kpath.special_points', inputdict)      # e.g., ["G","X","W","K","G"]
npts_per_seg = int(jmespath.search('kpath.npoints_segment', inputdict))  # e.g., 51

# Build a single continuous path as a list of fractional q-points (in units of reciprocal lattice)
band_path = [[special_points[sym] for sym in labels_flat]]  # NOTE: extra [] → list of paths
labels = [labels_flat]  # match the nesting: one label-list per path

# Let Phonopy create the dense q-path (npoints is per segment)
qpoints, connections = get_band_qpoints_and_path_connections(
    band_path,
    npoints=npts_per_seg
)

# Compute and save band structure
phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
phonon.write_yaml_band_structure(filename='./phonopy_band.yaml')
phonon.write_hdf5_band_structure(filename='./phonopy_band.h5')
'''

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
                    'path':  'generator.dest_dir',
                    'value': './',
                },
                {
                    'path':  'generator.steps',
                    'value': ['scf_qe', 'scf_pw2bgw_qe'],
                },
                {
                    'path': 'manager.steps',
                    'value': ['scf_qe', 'scf_pw2bgw_qe'],
                },
                {
                    'path': 'manager.plots',
                    'value': ['scf_qe', 'scf_pw2bgw_qe'],
                },
                {
                    'path':  'manager.dest_dir',
                    'value': './',
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
            './band.yaml',
            './phonopy_band.yaml',
            './phonopy_band.h5',
            './phonopy',
            './phonopy.yaml',
            './phonopy_disp.yaml',
        ]
    
    def plot(self):
        PhonopyPlot().save_figures()

#endregion
