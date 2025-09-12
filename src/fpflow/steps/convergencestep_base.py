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
class ConvergenceBaseStep(Step):
    def __init__(self, generatorclass=None, subdir1:str = None, subdir2: str = None, **kwargs):
        '''
        Problem:
            Generator calls Step which calls Generator. 
            To lazy to figure out and implement dependency inversion. So,
            I pass the class to the contructor and choose composition as a design choice. 
        '''
        super().__init__(**kwargs)
        self.generatorclass = generatorclass
        self.subdir1: str = subdir1
        self.subdir2: str = subdir2

        assert len(self.subdir1)>0, 'subdir1 should have length more than 1'
        assert len(self.subdir2)>0, 'subdir2 should have length more than 1'

    @property
    def script_conv(self) -> str:
        return f'''import os
import sys
import glob 
import subprocess
from fpflow.managers.run import subprocess_run

# Set the stdout and stderr. 
outfile = open('script_conv{self.subdir1}{self.subdir2}.py.out', 'w')
sys.stdout = outfile
sys.stderr = outfile

# Get the directories. 
{self.subdir1}{self.subdir2}_dirs = [inode for inode in glob.glob('./convergence/{self.subdir1}/{self.subdir2}/*') if os.path.isdir(inode)]
{self.subdir1}{self.subdir2}_dirs.sort()
start: int = 0
stop: int = len({self.subdir1}{self.subdir2}_dirs)

# Override if needed. Comment this out and set. 
#start = 0 
#stop = 1

total_time: float = 0.0
for {self.subdir1}{self.subdir2}_dir in {self.subdir1}{self.subdir2}_dirs[start:stop]:
    total_time = subprocess_run('./run.sh', total_time=total_time, dest_dir={self.subdir1}{self.subdir2}_dir)

print(f'Done conv{self.subdir1}{self.subdir2} in total time: ', total_time, ' seconds.', flush=True)
'''
    
    @property
    def job_conv(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, f'convergence.{self.subdir1}.{self.subdir2}.job_info')
        
        return f'''#!/bin/bash
{scheduler.get_script_header()}

./link.sh

python ./script_conv{self.subdir1}{self.subdir2}.py
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
        params: list = jmespath.search(f'convergence.{self.subdir1}.{self.subdir2}.dir_list[{dir_list_idx}].params', inputdict_local)
        link_inodes: list = jmespath.search(f'convergence.{self.subdir1}.{self.subdir2}.dir_list[{dir_list_idx}].link_inodes', inputdict_local)
        generator_steps: list = jmespath.search(f'convergence.{self.subdir1}.{self.subdir2}.dir_list[{dir_list_idx}].generator', inputdict_local)

        # Update values. 
        for param in params:
            glom.assign(inputdict_local, param['path'], param['value'])
        glom.assign(inputdict_local, 'generator.dest_dir', './')
        glom.assign(inputdict_local, 'generator.pre_steps', generator_steps['pre_steps'])
        glom.assign(inputdict_local, 'generator.steps', generator_steps['steps'])
        glom.assign(inputdict_local, 'generator.post_steps', generator_steps['post_steps'])
        glom.assign(inputdict_local, 'manager.steps', generator_steps['steps'])

        # Set link files. 
        if link_inodes is not None and len(link_inodes)>0:
            for link_inode in link_inodes:
                os.system(f'ln -sf {link_inode["source"]} {link_inode["dest"]}')

        # Write input yaml file. 
        InputYaml.to_yaml_file('./input.yaml', inputdict_local)

        # Run generator in the subdirectory. 
        generator = self.generatorclass.from_inputyaml('./input.yaml')
        generator.create()

    def create_in_subdirs(self):
        dir_list_len = len(jmespath.search(f'convergence.{self.subdir1}.{self.subdir2}.dir_list', self.inputdict))

        for dir_list_idx in range(dir_list_len):
            self.dest_dir: str = f'./convergence/{self.subdir1}/{self.subdir2}/dset_{dir_list_idx}'
            self.current_dir: str = os.getcwd()
            os.makedirs(self.dest_dir, exist_ok=True)
            self._create_in_subdir(dir_list_idx)
                
    def create(self):
        self.create_in_subdirs()

        extra_filecontents: dict = {
            f'script_conv{self.subdir1}{self.subdir2}.py': self.script_conv,
            f'job_conv{self.subdir1}{self.subdir2}.sh': self.job_conv,
        }
        for filename, filecontents in extra_filecontents.items():
            str_2_f(filecontents, filename)
        
        os.system('chmod u+x ./*.sh')

    @property
    def job_scripts(self) -> List[str]:
        return [
            f'./job_conv{self.subdir1}{self.subdir2}.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './convergence',
            f'./script_conv{self.subdir1}{self.subdir2}.py.out',
            f'./script_conv{self.subdir1}{self.subdir2}.py',
            f'./job_conv{self.subdir1}{self.subdir2}.sh',
        ]

#endregion