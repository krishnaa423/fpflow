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
class DmcXctphStep(Step):
    @property
    def dmc_xctph(self) -> str:
        filestring = f'''#!/bin/bash
from fpflow.analysis.dmc_xctph.dmc_xctph import DmcXctph
dmc_xctph: DmcXctph = DmcXctph()
dmc_xctph.run()
dmc_xctph.write()
'''
        return filestring

    @property
    def job_dmc_xctph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dmc_xctph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}python ./dmc_xctph.py &> dmc_xctph.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './dmc_xctph/dmc_xctph.py': self.dmc_xctph,
            './dmc_xctph/job_dmc_xctph.sh': self.job_dmc_xctph,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './dmc_xctph/job_dmc_xctph.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dmc_xctph',
        ]
    
#endregion