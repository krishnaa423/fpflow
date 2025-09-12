#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.nestedstep_base import NestedBaseStep
import jmespath 
from fpflow.schedulers.scheduler import Scheduler
import copy 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class MlQeDftStep(NestedBaseStep):
    def __init__(self, **kwargs):
        super().__init__(steptag='mldftqe', **kwargs)

    @property
    def folder_inputdict_changes_map(self) -> dict:
        '''
        Return:
            {
                'folder1': [
                    {
                        'path': 'some.path1',
                        'value': 12,
                    },
                    {
                        'path': 'some.path2',
                        'value': 45,
                    },
                ],
                'folder2': [
                    {
                        'path': 'some.path1',
                        'value': 123,
                    },
                    {
                        'path': 'some.path2',
                        'value': 345,
                    },
                ],
            }
        '''
        num_structures = len(jmespath.search('structures.list[*]', self.inputdict))

        output = {}

        # Exclude steps that nested steps. Or else it becomes infinitely nested.  
        generator_steps: list = [item for item in copy.deepcopy(jmespath.search('generator.steps', self.inputdict)) if not issubclass(self.stepmap[item], NestedBaseStep)]
        manager_steps: list = [item for item in copy.deepcopy(jmespath.search('manager.steps', self.inputdict)) if not issubclass(self.stepmap[item], NestedBaseStep)]

        for struct_idx in range(num_structures):
            output[f'./mldftqe/dset_{struct_idx}/'] = [
                {
                    'path': 'structures.active_idx',
                    'value': struct_idx,
                },
                {
                    'path': 'generator.dest_dir',
                    'value': './',
                },
                {
                    'path': 'generator.steps',
                    'value': generator_steps,
                },
                {
                    'path': 'manager.steps',
                    'value': manager_steps,
                },
            ]

        return output

    @property
    def extra_filecontents(self) -> dict:
        output: dict = super().extra_filecontents

        return output
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_mldftqe.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './mldftqe',
            './job_mldftqe.sh',
            './script_mldftqe.py',
            './script_mldftqe.py.out',
        ]

#endregion