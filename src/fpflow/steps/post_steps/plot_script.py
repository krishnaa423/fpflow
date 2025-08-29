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
class PlotScriptStep(Step):
    @property
    def file_contents(self) -> dict:
        return {
            'plot.sh': f'''#!/usr/bin/env python3

from fpflow.managers.manager import Manager

manager = Manager.from_inputyaml('./input.yaml')
manager.plot()
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
        return ['plot.sh']
#endregion