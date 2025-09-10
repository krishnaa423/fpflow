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
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePhmodesStep(Step):
    @property
    def dynmat(self) -> str:
        qpt_idx = jmespath.search('phmodes.qpt_idx', self.inputdict)
        dynmatdict: dict = {
            'input': {
                'asr': "'crystal'",
                'fildyn': f"'struct.dyn{qpt_idx}'",
                'filxsf': "'struct_phmodes.axsf'",
            }
        }

        # Update if needed. 
        update_dict(dynmatdict, jmespath.search('phmodes.args', self.inputdict))

        return NamelistGrammar().write(dynmatdict)
    
    @property
    def job_dynmat(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'phmodes.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}dynmat.x < dynmat.in &> dynmat.in.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'dynmat.in': self.dynmat,
            'job_dynmat.sh': self.job_dynmat,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_dynmat.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dynmat.in',
            './dynmat.out',
            './dynmat.in.out',
            './dynmat.mold',
            './input_tmp.in',
            './job_dynmat.sh',
            './struct.dyn*',
            './struct_phmodes.axsf',
        ]
#endregion