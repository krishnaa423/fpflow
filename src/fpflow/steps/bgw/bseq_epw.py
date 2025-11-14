#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
import copy 
import numpy as np 
import jmespath
from fpflow.structure.kpts import Kpts
from fpflow.steps.bgw.kernel import BgwKernelStep
from fpflow.steps.bgw.absorption import BgwAbsorptionStep
from fpflow.steps.bgw.plotxct import BgwPlotxctStep
from fpflow.schedulers.scheduler import Scheduler
import glom 
from fpflow.inputs.inputyaml import InputYaml
from fpflow.io.change_dir import change_dir
import glob
import h5py

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwBseqStep(Step):
    def __init__(self, generatorclass=None, **kwargs):
        '''
        Problem:
            Generator calls Step which calls Generator. 
            To lazy to figure out and implement dependency inversion. So,
            I pass the class to the contructor and choose composition as a design choice. 
        '''
        super().__init__(**kwargs)
        self.generatorclass = generatorclass

    def get_xct_h5(self):
        # Get sizes.
        nQ: int = len(glob.glob('./q_*'))
        with h5py.File('./q_0/eigenvectors.h5', 'r') as h5file:
            nS: int = h5file['/exciton_header/kpoints/nevecs'][()]
            nk: int = h5file['/exciton_header/kpoints/nk'][()]
            nc: int = h5file['/exciton_header/params/nc'][()]
            nv: int = h5file['/exciton_header/params/nv'][()]
            ns: int = h5file['/exciton_header/params/ns'][()]

        # Create dataset arrays.
        eigs: np.ndarray = np.zeros(shape=(nQ, nS), dtype='f8')
        evecs: np.ndarray = np.zeros(shape=(nQ, nS, nk, nc, nv), dtype='c16')

        for qpt_idx in range(nQ):
            with h5py.File(f'./q_{qpt_idx}/eigenvectors.h5', 'r') as h5file:
                eigs[qpt_idx, :] = h5file['/exciton_data/eigenvalues'][:]
                evecs[qpt_idx, :, :, :, :] = h5file['/exciton_data/eigenvectors'][0, :, :, :, :, 0, 0] + 1j * h5file['/exciton_data/eigenvectors'][0, :, :, :, :, 0, 1]

        # Write to h5 file.
        dsets: dict = {
            '/eigs': eigs,
            '/evecs': evecs,
        }
        with h5py.File('./xct.h5', 'w') as h5file:
            for dset_name, dset_data in dsets.items():
                h5file.create_dataset(f'{dset_name}', data=dset_data)

    @property
    def script_bseq_epw(self):
        return f'''import os
import sys
import glob
import subprocess
from fpflow.managers.run import subprocess_run

# Set the stdout and stderr. 
outfile = open('script_bseq_epw.py.out', 'w')
sys.stdout = outfile
sys.stderr = outfile

# Get the directories. 
bseq_dirs = [inode for inode in glob.glob('./eigv/*') if os.path.isdir(inode)]
bseq_dirs.sort()
start: int = 0
stop: int = len(bseq_dirs)

# Override if needed. Comment this out and set. 
#start = 0 
#stop = 1

total_time: float = 0.0
for bseq_dir in bseq_dirs[start:stop]:
    total_time = subprocess_run('./job_kernel.sh', total_time=total_time, dest_dir=bseq_dir)
    total_time = subprocess_run('./job_absorption.sh', total_time=total_time, dest_dir=bseq_dir)
    total_time = subprocess_run('./job_plotxct.sh', total_time=total_time, dest_dir=bseq_dir)
    
print(f'Done bseq in total time: ', total_time, ' seconds.', flush=True)
'''

    @property
    def job_bseq_epw(self):
        bseq_scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bseq.job_info')

        return f'''#!/bin/bash
{bseq_scheduler.get_script_header()}

python ./script_bseq_epw.py

# Write xct.h5 file.
bseq_pp="
from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.bgw.bseq_epw import EpwBseqStep
inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
bseq: EpwBseqStep = EpwBseqStep(inputdict=inputdict)
bseq.get_xct_h5()
"
python -c "$bseq_pp" &> bseq_pp.out 
'''

    @change_dir
    def _create_in_subdir(self, Qpt: list, Qpt_idx: int):
        '''
        This function:
        - Changes to subdirectory. 
        - Creates the input.yaml file there. 
        - Runs the generator to create the step files. 
        '''
        inputdict_local: dict = copy.deepcopy(self.inputdict)

        # Change parameters. 
        glom.assign(inputdict_local, 'bse.kernel.qshift', Qpt)
        glom.assign(inputdict_local, 'bse.kernel.link_dir_prefix', '../../')
        glom.assign(inputdict_local, 'bse.absorption.qshift', Qpt)
        glom.assign(inputdict_local, 'bse.absorption.link_dir_prefix', '../../')
        glom.assign(inputdict_local, 'bse.plotxct.link_dir_prefix', '../../')
        # Change generator and manager parameters. 
        glom.assign(inputdict_local, 'generator.dest_dir', './')
        glom.assign(inputdict_local, 'generator.pre_steps', [])
        glom.assign(inputdict_local, 'generator.steps', ['kernel_bgw', 'absorption_bgw', 'plotxct_bgw'])
        glom.assign(inputdict_local, 'generator.post_steps', [
            'create_script',
            'run_script',
            'remove_script',
            'plot_script',
            'interactive_script',
        ])
        # glom.assign(inputdict_local, 'manager.steps', ['kernel_bgw', 'absorption_bgw', 'plotxct_bgw'])
        glom.assign(inputdict_local, 'manager.steps', ['kernel_bgw', 'absorption_bgw']) # Skip plotxct for now.

        # Write input yaml file. 
        InputYaml.to_yaml_file('./input.yaml', inputdict_local)

        # Run generator in the subdirectory. 
        generator = self.generatorclass.from_inputyaml('./input.yaml')
        generator.create()

    def create_in_subdirs(self):
        Qgrid: np.ndarray = jmespath.search('bseq.qgrid', self.inputdict)
        is_reduced: bool = jmespath.search('bseq.sym', self.inputdict)
        Qpts: list = Kpts.from_kgrid(
            kgrid=Qgrid,
            is_reduced=is_reduced,
        ).bseq_epw_qpts

        for Qpt_idx, Qpt in enumerate(Qpts):
            self.dest_dir: str = f'./eigv/q_{Qpt_idx}'
            self.current_dir: str = os.getcwd()
            os.makedirs(self.dest_dir, exist_ok=True)
            self._create_in_subdir(Qpt, Qpt_idx)

    def create(self):
        self.create_in_subdirs()

        extra_filecontents: dict = {
            'script_bseq_epw.py': self.script_bseq_epw,
            'job_bseq_epw.sh': self.job_bseq_epw,
        }

        for filename, filecontents in extra_filecontents.items():
            str_2_f(filecontents, filename)
        
        os.system('chmod u+x ./*.sh')

    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_bseq_epw.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []

    @property
    def remove_inodes(self) -> List[str]:
        return [
            './eigv',
            './script_bseq_epw.py',
            './script_bseq_epw.py.out',
            './job_bseq_epw.sh',
        ]

#endregion