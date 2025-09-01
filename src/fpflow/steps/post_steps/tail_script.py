#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class TailScriptStep(Step):
    @property
    def file_contents(self) -> dict:
        return {
            'tail.sh': f'''#!/bin/bash
tail -n100 -f $1
'''
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return []

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return ['./follow']
#endregion