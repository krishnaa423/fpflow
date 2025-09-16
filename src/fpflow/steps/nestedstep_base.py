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
from fpflow.inputs.inputyaml import InputYaml
from fpflow.io.change_dir import change_dir
from abc import ABC, abstractmethod
from benedict import benedict
import yaml 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class NestedBaseStep(Step, ABC):
    def __init__(self, generatorclass=None, steptag:str=None, stepmap:dict=None, **kwargs):
        '''
        Problem:
            Generator calls Step which calls Generator. 
            To lazy to figure out and implement dependency inversion. So,
            I pass the class to the contructor and choose composition as a design choice. 
        '''
        super().__init__(**kwargs)
        self.generatorclass = generatorclass
        self.steptag: str = steptag
        self.stepmap: dict = stepmap

    @property
    def script_steptag(self):
        return f'''import os
import sys
import glob
import subprocess
from fpflow.managers.run import subprocess_run

# # Set the stdout and stderr. 
# outfile = open('script_{self.steptag}.py.out', 'w')
# sys.stdout = outfile
# sys.stderr = outfile

# Get the directories. 
{self.steptag}_dirs = [inode for inode in glob.glob('./{self.steptag}/*') if os.path.isdir(inode)]
{self.steptag}_dirs.sort()
start: int = 0
stop: int = len({self.steptag}_dirs)

# Override if needed. Comment this out and set. 
#start = 0 
#stop = 1

total_time: float = 0.0
for {self.steptag}_dir in {self.steptag}_dirs[start:stop]:
    total_time = subprocess_run('./run.sh', total_time=total_time, dest_dir={self.steptag}_dir)

print(f'Done {self.steptag} in total time: ', total_time, ' seconds.', flush=True)
'''

    @property
    def job_steptag(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, f'{self.steptag}.job_info')

        return f'''#!/bin/bash
{scheduler.get_script_header()}

python ./script_{self.steptag}.py &> script_{self.steptag}.py.out
'''

    @property
    @abstractmethod
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
        return {}
    
    @property
    @abstractmethod
    def extra_filecontents(self) -> dict:
        '''
        Return:
            {
                'filename1': contents1,
                'filename2': contents2,
            }
        '''
        return {
            f'job_{self.steptag}.sh': self.job_steptag,
            f'script_{self.steptag}.py': self.script_steptag,
        }

    @change_dir
    def _create_in_subdir(self, dirname: str):
        '''
        This function:
        - Changes to subdirectory. 
        - Creates the input.yaml file there. 
        - Runs the generator to create the step files. 
        '''
        inputdict_local = benedict(copy.deepcopy(self.inputdict))
        input_changes_list: list = self.folder_inputdict_changes_map[dirname]

        # Update values. 
        for change in input_changes_list:
            inputdict_local.set(change['path'], change['value'])

        # Write input yaml file. 
        yaml.safe_dump(dict(inputdict_local), open("input.yaml", "w"))

        # Run generator in the subdirectory. 
        generator = self.generatorclass.from_inputyaml('./input.yaml')
        generator.create()

    def create_in_subdirs(self):
        dirnames: list = self.folder_inputdict_changes_map.keys()

        for dirname in dirnames:
            self.dest_dir: str = dirname
            self.current_dir: str = os.getcwd()
            os.makedirs(self.dest_dir, exist_ok=True)
            self._create_in_subdir(dirname)
                
    def create(self):
        # Create step files in sub directories.
        self.create_in_subdirs()

        for filename, filecontents in self.extra_filecontents.items():
            # Write content to file. 
            str_2_f(filecontents, filename)

            # Change permissions if executable script. 
            if filename.endswith('.sh'):
                os.system(f'chmod u+x {filename}')

    @property
    @abstractmethod
    def job_scripts(self) -> List[str]:
        return []

    @property
    @abstractmethod
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    @abstractmethod
    def remove_inodes(self) -> List[str]:
        return []

#endregion