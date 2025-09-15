#region modules
from typing import List
from fpflow.steps.step import Step 
from pathlib import Path
from fpflow.io.change_dir import change_dir
from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.steps_map import step_class_map
import os 
import jmespath
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class Manager:
    def __init__(self, **kwargs):
        self.inputdict: dict = None 
        self.steps: List[Step] = None 
        self.plot_steps: List[Step] = None 
        self.current_dir: str = None 
        self.dest_dir: str = None 
        self.is_recursive: bool = False 
        self.total_time: float = 0.0

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

        # Get data. 
        steps_list_str: List[str] = jmespath.search('manager.steps[*]', inputdict) if jmespath.search('manager.steps[*]', inputdict) is not None else []
        plot_steps_list_str: List[str] = jmespath.search('manager.plots[*]', inputdict) if jmespath.search('manager.plots[*]', inputdict) is not None else []
        steps: List[Step] = cls._get_step_classes(inputdict, steps_list_str)
        plot_steps: List[Step] = cls._get_step_classes(inputdict, plot_steps_list_str)
        dest_dir: str = jmespath.search('manager.dest_dir', inputdict)
        is_recursive: str = jmespath.search('manager.is_recursive', inputdict)

        return cls(
            inputdict=inputdict, 
            current_dir=os.getcwd(),
            dest_dir=dest_dir,
            steps=steps,
            plot_steps=plot_steps,
            is_recursive=is_recursive,
        )

    @change_dir
    def plot(self, **kwargs):
        for step in self.plot_steps:
            step.plot(**kwargs)

    @change_dir
    def run(self):
        self.total_time = 0.0
        for step in self.steps:
            self.total_time = step.run(total_time=self.total_time)

        return self.total_time

#endregion