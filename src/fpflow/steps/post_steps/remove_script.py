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
class RemoveScriptStep(Step):
    @property
    def file_contents(self) -> dict:
        return {
            'remove.sh': f'''#!/usr/bin/env python3

from fpflow.generators.generator import Generator

generator = Generator.from_inputyaml('./input.yaml')
generator.remove()
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
        return ['./remove.sh']
#endregion