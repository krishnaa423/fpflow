#region modules
from argparse import ArgumentParser 
from importlib.util import find_spec
import os 
from fpflow.generators.generator import Generator
#endregion

#region variables
#endregion

#region functions
def fpflow():
    parser = ArgumentParser(description="fpflow templates, generator, and manager.")
    parser.add_argument('--input', type=str, nargs='?', const='input_Si.yaml', default=None, help='Provide a template to genrate input yaml file form.')
    parser.add_argument('--generator', nargs='?', const='create', default=None, type=str, help='Create or remove files')
    parser.add_argument('--manager', nargs='?', const='interactive', default=None, type=str, help='Run the run.sh script or rund.sh for interactive or background manager runs.')
    args = parser.parse_args()

    if args.input!=None:
        pkg_dir = os.path.dirname(find_spec('fpflow').origin)
        filename = os.path.join(pkg_dir, 'data', args.input)
        os.system(f'cp {filename} ./input.yaml')

    if args.generator!=None:
        match args.generator:
            case 'create':
                generator = Generator.from_inputyaml()
                generator.create()
            case 'remove':
                generator = Generator.from_inputyaml()
                generator.remove()
            case _:
                raise ValueError('generator needs to have create or remove as arguments')

    if args.manager!=None:
        match args.run_manager:
            case 'interactive':
                os.system('./run.sh')
            case 'background':
                os.system('./rund.sh')
            case _:
                raise ValueError('run_manager needs to have interactive or background as arguments')
            
#endregion

#region classes
#endregion