#region modules
from importlib.util import find_spec
import os 
import yaml 
import jmespath
from fpflow.schedulers.jobinfo import JobInfo
from fpflow.io.logging import get_logger
from fpflow.inputs.inputyaml import InputYaml
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
    def from_scheds_dict(cls, system_name: str, scheds_dict: dict):
        system_dict: dict = jmespath.search(f'{system_name}', scheds_dict)
        is_interactive: bool = False if jmespath.search('is_interactive', system_dict) is None else jmespath.search('is_interactive', system_dict)

        match is_interactive:
            case False:
                return Scheduler(**system_dict)
            case True:
                return Interactive(**system_dict)

    @classmethod
    def from_sched_yamlfilename(cls, system_name: str, filename: str = None, update_dict: dict={}):
        if filename is None:
            pkg_dir = os.path.dirname(find_spec('fpflow').origin)
            filename = os.path.join(pkg_dir, 'data', 'schedulers.yaml')

        scheds_dict: dict = None 
        with open(filename, 'r') as f: scheds_dict = yaml.safe_load(f)
        
        system_dict: dict = jmespath.search(f'{system_name}', scheds_dict)
        if system_dict is None: system_dict = {}
        system_dict.update(update_dict)
        scheds_dict.update({system_name: system_dict})

        return cls.from_scheds_dict(system_name, scheds_dict)
    
    @classmethod
    def from_inputdict(cls, inputdict: dict):
        system_name: str = jmespath.search('scheduler.manager', inputdict)
        update_dict: dict = {} if jmespath.search('scheduler.args', inputdict) is None else jmespath.search('scheduler.args', inputdict)

        return cls.from_sched_yamlfilename(system_name=system_name, update_dict=update_dict)

    @classmethod
    def from_input_yamlfilename(cls, filename: str='./input.yaml'):
        inputdict: dict = InputYaml.from_yaml_file(filename).inputdict

        return cls.from_inputdict(inputdict)
        
    def get_script_header(self, info: JobInfo):
        header =  ''

        for key, value in self.core_header.items():
            header += f'#{self.batchscript_exec.upper()} --{key}={value}\n'

        if len(self.core_header.keys())>0:
            header += f'#{self.batchscript_exec.upper()} --nodes={info.nodes}\n'
            header += f'#{self.batchscript_exec.upper()} --time={info.time}\n'

        header += f'{"" if self.extra_commands is None else self.extra_commands}\n'

        return header 

    def get_exec_prefix(self, info: JobInfo):
        prefix = ''

        if self.mpi_exec is not None and self.mpi_exec!='':
            prefix += f'{self.mpi_exec} -n {info.ntasks} '

            if self.is_gpu:
                prefix += f'--gpus-per-task={jmespath.search("gpus", self.node_info)} '

        return prefix 
        
    def get_exec_infix(self, info: JobInfo):
        infix = ''

        if self.mpi_exec is not None and self.mpi_exec!='':
            if info.nk is not None:
                infix += f'-nk {info.nk} '
            
            if info.ni is not None:
                infix += f'-ni {info.ni} '

        return infix 

    def get_launch_exec(self, info: JobInfo):
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

    def get_script_header(self, info: JobInfo):
        header =  ''

        header += f'{"" if self.extra_commands is None else self.extra_commands}\n'

        return header 

    def get_exec_prefix(self, info: JobInfo):
        return super().get_exec_prefix(info)
        
    def get_exec_infix(self, info: JobInfo):
        return super().get_exec_infix(info)

    def get_launch_exec(self, info: JobInfo):
        return super().get_launch_exec(info)
    
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