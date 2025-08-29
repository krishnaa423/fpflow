#region modules
from fpflow.steps.step import Step
from fpflow.steps.steps_map import step_class_map
from typing import List 
from fpflow.inputs.inputyaml import InputYaml
import jmespath
import os
import functools
from fpflow.io.change_dir import change_dir
from fpflow.io.logging import get_logger
#endregion

#region variables
logger = get_logger()
#endregion

#region functions

#endregion

#region classes
class Generator:
    '''
    Typical usage:
        generator = Generator.from_inputyaml('./input.yaml')
        generator.doall()
    '''
    def __init__(self, **kwargs):
        self.inputdict: dict = None 
        self.presteps: List[Step] = None 
        self.steps: List[Step] = None 
        self.poststeps: List[Step] = None 
        self.current_dir: str = None 
        self.dest_dir: str = None 

        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def _get_step_classes(cls, inputdict: dict, steps_list_str: List[str]) -> List[Step]:
        step_classes = []
        for step_str in steps_list_str:
            if step_str in step_class_map.keys():
                step_classes.append(step_class_map[step_str](inputdict=inputdict))
            else:
                print(f'{step_str} does not have a class map value.', flush=True)
        
        return step_classes

    @classmethod
    def from_inputyaml(cls, filename: str='./input.yaml'):
        # Copy first and create a filename that is standard. 
        if filename!='./input.yaml':
            os.system(f'cp {filename} ./input.yaml')

        # Create input dict. 
        inputdict: dict = InputYaml.from_yaml_file(filename).inputdict

        # Get the steps from yaml. 
        presteps_list_str: List[str] = jmespath.search('generator.pre_steps[*]', inputdict) if jmespath.search('generator.pre_steps[*]', inputdict)!=None else []
        steps_list_str: List[str] = jmespath.search('generator.steps[*]', inputdict) if jmespath.search('generator.steps[*]', inputdict)!=None else []
        poststeps_list_str: List[str] = jmespath.search('generator.post_steps[*]', inputdict) if jmespath.search('generator.post_steps[*]', inputdict)!=None else []

        # Get the step classes. 
        presteps: List[Step] = cls._get_step_classes(inputdict, presteps_list_str) 
        steps: List[Step] = cls._get_step_classes(inputdict, steps_list_str) 
        poststeps: List[Step] = cls._get_step_classes(inputdict, poststeps_list_str)

        dest_dir: str = jmespath.search('generator.dest_dir', inputdict)

        return cls(
            inputdict=inputdict, 
            current_dir=os.getcwd(),
            dest_dir=dest_dir,
            presteps=presteps, 
            steps=steps,
            poststeps=poststeps,
        )
    
    @change_dir
    def create_presteps(self):
        for prestep in self.presteps:
            prestep.create()

    @change_dir
    def run_presteps(self):
        for prestep in self.presteps:
            prestep.run()

    @change_dir
    def create_steps(self):
        for step in self.steps:
            step.create()

    @change_dir
    def create_poststeps(self):
        for poststep in self.poststeps:
            poststep.create()

    @change_dir
    def run_poststeps(self):
        for poststep in self.poststeps:
            poststep.run()

    def create(self):
        self.create_presteps()
        self.run_presteps()
        self.create_steps()
        self.create_poststeps()
        self.run_poststeps()

    @change_dir
    def remove_presteps(self):
        for prestep in self.presteps:
            prestep.remove()

    @change_dir
    def remove_steps(self):
        for step in self.steps:
            step.remove()

    def remove(self):
        self.remove_presteps()
        self.remove_steps()

#endregion