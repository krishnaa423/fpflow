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
class EsrgenStep(Step):
    @property
    def file_contents(self) -> dict:
        return {}
    
    @property
    def job_scripts(self) -> List[str]:
        return []

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return []
    
    def __init__(self, inputdict: dict, **kwargs):
        self.inputdict: dict = inputdict 
#endregion