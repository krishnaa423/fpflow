#region modules
import yaml  
from fpflow.io.read_write import str_2_f
import jmespath 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class InputYaml:
     def __init__(self, inputdict: dict):
         self.inputdict: dict = inputdict  

     @classmethod
     def from_yaml_file(cls, yaml_filename: str='./input.yaml'):
          with open(yaml_filename, 'r') as f: inputdict: dict = yaml.safe_load(f)
          return InputYaml(inputdict)

     @classmethod
     def from_yaml_str(cls, yaml_string: str):
          inputdict: dict = yaml.safe_load(yaml_string)
          return InputYaml(inputdict)

     @classmethod
     def to_yaml_str(cls, inputdict: dict) -> str:
          filestring: str = yaml.dump(inputdict)
          return filestring

     @classmethod 
     def to_yaml_file(cls, filename: str, inputdict: dict):
          str_2_f(cls.to_yaml_str(inputdict), filename)
    
#endregion