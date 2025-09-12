#region modules
import yaml 
import jmespath
from glom import assign
import os 
#endregion

#region variables
#endregion

#region functions
def str_2_f(file_contents: str, filename: str):
    if filename is not None and os.path.dirname(filename) != '':
        os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as f: f.write(file_contents)

def f_2_str(filename: str) -> str:
    with open(filename, 'r') as f: return f.read()

#endregion

#region classes
#endregion