#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.schedulers.scheduler import Scheduler

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePpElphEpwStep(Step):
    @property
    def elph(self) -> str:
        file_string: str = f'''
#region modules
from fpelph.epw import Epw
#endregion

#region variables
#endregion

#region functions
def main():
    elph: Epw = Epw()
    elph.read()
    elph.write()
#endregion

#region classes
#endregion

#region main
main()
#endregion
'''
        return file_string
    
    @property
    def job_elph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'elph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf elph.out
touch elph.out
exec &> elph.out

{scheduler.get_exec_prefix()}python script_elph.py &> script_elph.py.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'script_elph.py': self.elph,
            'job_elph.sh': self.job_elph,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_elph.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './script_elph.py',
            './script_elph.py.out',
            './job_elph.sh',
            './elph.h5',
        ]
    
#endregion