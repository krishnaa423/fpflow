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
class SteStep(Step):
    @property
    def ste(self) -> str:
        filestring = f'''#!/bin/bash
from fpflow.analysis.ste.ste import Ste
ste: Ste = Ste()
ste.run()
ste.write()
'''
        return filestring

    @property
    def job_ste(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'ste.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}python ./ste.py &> ste.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './ste/ste.py': self.ste,
            './ste/job_ste.sh': self.job_ste,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './ste/job_ste.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './ste',
        ]
    
#endregion