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
class XctpolStep(Step):
    @property
    def xctpol(self) -> str:
        filestring = f'''#!/bin/bash
from fpflow.analysis.xctpol.xctpol import Xctpol
xctpol: Xctpol = Xctpol()
xctpol.run()
xctpol.write()
'''
        return filestring

    @property
    def job_xctpol(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'xctpol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}python ./xctpol.py &> xctpol.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './xctpol/xctpol.py': self.xctpol,
            './xctpol/job_xctpol.sh': self.job_xctpol,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './xctpol/job_xctpol.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './xctpol',
        ]
    
#endregion