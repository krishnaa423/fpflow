#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.managers.run import subprocess_run
import glob 
from abc import ABC, abstractmethod
#endregion

#region variables

#endregion

#region functions
#endregion

#region classes
class Step:
    @property
    def file_contents(self) -> dict:
        return {}
    
    @property
    def job_scripts(self) -> List[str]:
        return []

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return []
    
    def __init__(self, inputdict: dict, **kwargs):
        self.inputdict: dict = inputdict 

    def create(self):
        for filename, filecontents in self.file_contents.items():
            str_2_f(filecontents, filename)

            # Set permissions for the destination file. 
            if filename[-3:]=='.sh':
                os.system(f'chmod u+x {filename}')

    def run(self, **kwargs):
        total_time: float = 0.0 if not 'total_time' in kwargs.keys() else kwargs['total_time']

        for script in self.job_scripts:
            total_time = subprocess_run(script, **kwargs)

        return total_time

    def save(self, foldername: str):
        for inode in self.save_inodes:
            os.system(f'mv {inode} {foldername}/')

    def remove(self):
        for inode in self.remove_inodes:
            os.system(f'rm -rf {inode}')

    def plot(self, **kwargs):
        pass
    
#endregion