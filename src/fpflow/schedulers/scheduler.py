#region modules
from importlib.util import find_spec
import os 
import yaml 
import jmespath
from fpflow.io.logging import get_logger
from fpflow.inputs.inputyaml import InputYaml
import importlib
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class Scheduler:

    def __init__(self, **kwargs):
        self.node_info: dict = None 
        self.header: dict = None 
        self.launch_info: dict = None 
        self.is_gpu: bool = None 
        self.is_interactive: bool = None 
        self.nodes: int = None 
        self.ntasks: int = None
        self.time: str = None 
        self.nk: int = None 
        self.ni: int = None 

        for key, value in kwargs.items():
            setattr(self, key, value)

        self.set_additional()

    def set_additional(self):
        self.gpus: int = jmespath.search('gpus', self.node_info)
        self.core_header: dict = {k: v for k, v in jmespath.search('non_interactive', self.header).items() if k!='extra_commands'}
        self.extra_commands: str = jmespath.search('non_interactive.extra_commands', self.header)
        self.mpi_exec: str = jmespath.search('non_interactive.mpi_exec', self.launch_info)
        self.batchscript_exec: str = jmespath.search('non_interactive.batchscript_exec', self.launch_info)
    
    @classmethod
    def from_jmespath(cls, inputdict: dict, jmspath: str):
        # Read scheduler info. 
        sched_dict: dict = {} if jmespath.search(jmspath, inputdict) is None else jmespath.search(jmspath, inputdict)

        # Update master scheduler info. 
        pkg_dir = os.path.dirname(find_spec('fpflow').origin)
        filename = os.path.join(pkg_dir, 'data', 'schedulers.yaml')
        schedulers_dict: dict = None 
        with open(filename, 'r') as f: schedulers_dict  = yaml.safe_load(f)

        schedulers_dict[jmespath.search('name', sched_dict)].update(sched_dict)

        # Get master scheduler info. 
        master_info = jmespath.search(jmespath.search('name', sched_dict), schedulers_dict)
        is_interactive = False if jmespath.search('is_interactive', master_info) is None else jmespath.search('is_interactive', master_info)

        # Call the construtor.
        if is_interactive is not None and is_interactive:
            return Interactive(**master_info)
        else:
            return Scheduler(**master_info)

        
    def get_script_header(self):
        header =  ''

        for key, value in self.core_header.items():
            header += f'#{self.batchscript_exec.upper()} --{key}={value}\n'

        if len(self.core_header.keys())>0:
            header += f'#{self.batchscript_exec.upper()} --nodes={self.nodes}\n'
            header += f'#{self.batchscript_exec.upper()} --time={self.time}\n'

        header += f'{"" if self.extra_commands is None else self.extra_commands}\n'

        return header 

    def get_exec_prefix(self):
        prefix = ''

        if self.mpi_exec is not None and self.mpi_exec!='':
            prefix += f'{self.mpi_exec} -n {self.ntasks} '

            if self.is_gpu:
                prefix += f'--gpus-per-task={jmespath.search("gpus", self.node_info)} '

        return prefix 
        
    def get_exec_infix(self):
        infix = ''

        if self.mpi_exec is not None and self.mpi_exec!='':
            if self.nk is not None:
                infix += f'-nk {self.nk} '
            
            if self.ni is not None:
                infix += f'-ni {self.ni} '

        return infix 

    def get_launch_exec(self):
        submit = ''

        if self.batchscript_exec is not None and self.batchscript_exec!='':
            submit += f'{self.batchscript_exec} '

        return submit

class Interactive(Scheduler):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def set_additional(self):
        self.gpus: int = jmespath.search('gpus', self.node_info)
        self.core_header: dict = {k: v for k, v in jmespath.search('interactive', self.header).items() if k!='extra_commands'}
        self.extra_commands: str = jmespath.search('interactive.extra_commands', self.header)
        self.mpi_exec: str = jmespath.search('interactive.mpi_exec', self.launch_info)
        self.batchscript_exec: str = ''
        self.interactive_exec: str = jmespath.search('interactive.interactive_exec', self.launch_info)

    def get_script_header(self):
        header =  ''

        header += f'{"" if self.extra_commands is None else self.extra_commands}\n'

        return header 

    def get_exec_prefix(self):
        return super().get_exec_prefix()
        
    def get_exec_infix(self):
        return super().get_exec_infix()

    def get_launch_exec(self):
        return super().get_launch_exec()
    
    def get_interactive_script_str(self) -> str:
        option_string = ''

        if self.core_header is not None:
            for key, value in self.core_header.items():
                option_string += f' --{key}={value}'

        file_contents = f'''#!/bin/bash
{jmespath.search('interactive.interactive_exec', self.launch_info)} {option_string}
'''
        
        return file_contents

#endregion