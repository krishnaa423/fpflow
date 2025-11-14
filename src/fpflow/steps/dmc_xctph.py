#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
import jmespath
from fpflow.plots.dmc_xctph import DmcXctphPlot
from fpflow.schedulers.scheduler import Scheduler

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class DmcXctphStep(Step):
    def calc_dmc(self):
        '''
        Being lazy and writing dmc logic here for now.
        Have to move to analysis folder or seperate folder as the code grows. 
        '''
        pass

    @property
    def job_dmc_xctph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'xctph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

dmc_xctph="
from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.dmc_xctph import DmcXctphStep
inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
dmc_xctph: DmcXctphStep = DmcXctphStep(inputdict=inputdict)
dmc_xctph.calc_dmc()
"

python -c "$dmc_xctph" &> dmc_xctph.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './dmc_xctph/job_dmc_xctph.sh': self.job_dmc_xctph
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './dmc_xctph/job_dmc_xctph.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dmc_xctph'
        ]
    
    def plot(self, **kwargs):
        DmcXctphPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion