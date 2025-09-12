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
import copy 
import glom 
from fpflow.inputs.inputyaml import InputYaml
from fpflow.io.change_dir import change_dir
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeConvergenceScfStep(Step):
    def __init__(self, generatorclass=None, **kwargs):
        '''
        Problem:
            Generator calls Step which calls Generator. 
            To lazy to figure out and implement dependency inversion. So,
            I pass the class to the contructor and choose composition as a design choice. 
        '''
        super().__init__(**kwargs)
        self.generatorclass = generatorclass

    @property
    def script_convqescf(self) -> str:
        return f'''import os
import sys
import glob 
import subprocess
from fpflow.managers.run import subprocess_run

# Set the stdout and stderr. 
outfile = open('script_convqescf.py.out', 'w')
sys.stdout = outfile
sys.stderr = outfile

# Get the directories. 
qescf_dirs = [inode for inode in glob.glob('./convergence/qe/scf/*') if os.path.isdir(inode)]
start: int = 0
stop: int = len(qescf_dirs)

# Override if needed. Comment this out and set. 
#start = 0 
#stop = 1

total_time: float = 0.0
for qescf_dir in qescf_dirs[start:stop]:
    total_time = subprocess_run('./run.sh', total_time=total_time, dest_dir=qescf_dir)

print(f'Done convqescf in total time: ', total_time, ' seconds.', flush=True)
'''
    
    @property
    def job_convqescf(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'convergence.qe.scf.job_info')
        
        return f'''#!/bin/bash
{scheduler.get_script_header()}

./link.sh

python ./script_convqescf.py
'''

    @change_dir
    def _create_in_subdir(self, dir_list_idx: int):
        '''
        This function:
        - Changes to subdirectory. 
        - Creates the input.yaml file there. 
        - Runs the generator to create the step files. 
        '''
        inputdict_local: dict = copy.deepcopy(self.inputdict)
        params: list = jmespath.search(f'convergence.qe.scf.dir_list[{dir_list_idx}].params', inputdict_local)
        link_inodes: list = jmespath.search(f'convergence.qe.scf.dir_list[{dir_list_idx}].link_inodes', inputdict_local)
        common_link_inodes: list = jmespath.search(f'convergence.qe.scf.link_inodes', inputdict_local)

        # Update values. 
        for param in params:
            glom.assign(inputdict_local, param['path'], param['value'])
        glom.assign(inputdict_local, 'generator.dest_dir', './')
        glom.assign(inputdict_local, 'generator.pre_steps', ['pseudos_qe'])
        glom.assign(inputdict_local, 'generator.steps', ['scf_qe', 'dftelbands_qe'])
        glom.assign(inputdict_local, 'generator.post_steps', [
            'create_script',
            'run_script',
            'remove_script',
            'plot_script',
            'interactive_script',
        ])
        glom.assign(inputdict_local, 'manager.steps', ['scf_qe', 'dftelbands_qe'])

        # Write input yaml file. 
        InputYaml.to_yaml_file('./input.yaml', inputdict_local)

        # Run generator in the subdirectory. 
        generator = self.generatorclass.from_inputyaml('./input.yaml')
        generator.create()

    def create_in_subdirs(self):
        dir_list_len = len(jmespath.search('convergence.qe.scf.dir_list', self.inputdict))

        for dir_list_idx in range(dir_list_len):
            self.dest_dir: str = f'./convergence/qe/scf/dset_{dir_list_idx}'
            self.current_dir: str = os.getcwd()
            os.makedirs(self.dest_dir, exist_ok=True)
            self._create_in_subdir(dir_list_idx)
                
    def create(self):
        self.create_in_subdirs()

        extra_filecontents: dict = {
            'script_convqescf.py': self.script_convqescf,
            'job_convqescf.sh': self.job_convqescf,
        }
        for filename, filecontents in extra_filecontents.items():
            str_2_f(filecontents, filename)
        
        os.system('chmod u+x ./*.sh')

    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_convqescf.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './convergence',
            './script_convqescf.py.out',
            './script_convqescf.py',
            './job_convqescf.sh',
        ]
#endregion