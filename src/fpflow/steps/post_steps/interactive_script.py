#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.schedulers.scheduler import Scheduler
from fpflow.steps.step import Step 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class InteractiveScriptStep(Step):
    @property
    def file_contents(self) -> dict:
        return {
            'interactive.sh': Scheduler.from_jmespath(inputdict=self.inputdict, jmspath='job_types.interactive').get_interactive_script_str()
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return []

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return ['./interactive.sh']
#endregion