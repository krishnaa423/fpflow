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
import copy 
import glom 
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
        generator.create()
        generator.remove()
    '''
    @property
    def inputdict(self):
        return self._inputdict
    
    @inputdict.setter
    def inputdict(self, value):
        self._inputdict = value 
        
        # Run extract attributes to update related info. 
        self._extract_attributes()

    def __init__(self, inputdict: dict):
        self._inputdict: dict = inputdict 
        self._presteps: List[Step] = None 
        self._steps: List[Step] = None 
        self._poststeps: List[Step] = None 
        self.current_dir: str = None 
        self.dest_dir: str = None 

        # Update related info. 
        self._extract_attributes()

    def _extract_attributes(self):
        # Get the steps from yaml. 
        presteps_list_str: List[str] = jmespath.search('generator.pre_steps[*]', self._inputdict) if jmespath.search('generator.pre_steps[*]', self._inputdict) is not None else []
        steps_list_str: List[str] = jmespath.search('generator.steps[*]', self._inputdict) if jmespath.search('generator.steps[*]', self._inputdict) is not None else []
        poststeps_list_str: List[str] = jmespath.search('generator.post_steps[*]', self._inputdict) if jmespath.search('generator.post_steps[*]', self._inputdict) is not None else []

        # Get the step classes. 
        self._presteps: List[Step] = self._get_step_classes(self._inputdict, presteps_list_str) 
        self._steps: List[Step] = self._get_step_classes(self._inputdict, steps_list_str) 
        self._poststeps: List[Step] = self._get_step_classes(self._inputdict, poststeps_list_str)
        self.current_dir: str = os.getcwd()
        self.dest_dir: str = jmespath.search('generator.dest_dir', self._inputdict)

    def _get_step_classes(cls, inputdict: dict, steps_list_str: List[str]) -> List[Step]:
        step_classes = []
        for step_str in steps_list_str:
            if step_str in step_class_map.keys():
                step_classes.append(step_class_map[step_str](inputdict=inputdict, generatorclass=Generator, stepmap=step_class_map))
            else:
                print(f'{step_str} does not have a class map value.', flush=True)
        
        return step_classes

    @classmethod
    def from_inputyaml(cls, filename: str='./input.yaml'):
        # Create input dict. 
        inputdict: dict = InputYaml.from_yaml_file(filename).inputdict
        return cls(inputdict=inputdict)
    
    @change_dir
    def create_presteps(self):
        for prestep in self._presteps:
            # print(f'Creating {prestep} in dir: {os.getcwd()}', flush=True)
            prestep.create()

    @change_dir
    def run_presteps(self):
        for prestep in self._presteps:
            prestep.run()

    @change_dir
    def create_steps(self):
        for step in self._steps:
            # print(f'Creating {step} in dir: {os.getcwd()}', flush=True)
            step.create()

    @change_dir
    def create_poststeps(self):
        for poststep in self._poststeps:
            # print(f'Creating {poststep} in dir: {os.getcwd()}', flush=True)
            poststep.create()

    @change_dir
    def run_poststeps(self):
        for poststep in self._poststeps:
            poststep.run()

    def write_inputyaml_with_destdir_changed(self):
        if self.dest_dir is not None and os.path.abspath(self.dest_dir)!=os.getcwd():
            os.makedirs(self.dest_dir, exist_ok=True)
            updated_dict: dict = copy.deepcopy(self.inputdict)
            glom.assign(updated_dict, 'generator.dest_dir', './')
            InputYaml.to_yaml_file(os.path.join(self.dest_dir, 'input.yaml'), updated_dict)

    def create(self):
        self.write_inputyaml_with_destdir_changed()
        self.create_presteps()
        self.run_presteps()
        self.create_steps()
        self.create_poststeps()
        self.run_poststeps()

    @change_dir
    def remove_presteps(self):
        for prestep in self._presteps:
            prestep.remove()

    @change_dir
    def remove_steps(self):
        for step in self._steps:
            step.remove()

    def remove(self):
        self.remove_presteps()
        self.remove_steps()

        # Delete the destination directory only if it is not the current directory. 
        if os.path.abspath(self.dest_dir) != os.getcwd():
            os.system(f'rm -rf {self.dest_dir}')

#endregion