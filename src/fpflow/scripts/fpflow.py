#region modules
from argparse import ArgumentParser 
from importlib.util import find_spec
import os 
from fpflow.generators.generator import Generator
from fpflow.managers.manager import Manager
#endregion

#region variables
#endregion

#region functions
def run_input_template(args):
    pkg_dir = os.path.dirname(find_spec('fpflow').origin)
    filename = os.path.join(pkg_dir, 'data', 'templates', args.template)
    os.system(f'cp {filename} ./input.yaml')

def run_input_list(args):
    pkg_dir = os.path.dirname(find_spec('fpflow').origin)
    template_dir = os.path.join(pkg_dir, 'data', 'templates')
    for template in os.listdir(template_dir): print(template, flush=True)
    

def run_generator_create(args):
    generator = Generator.from_inputyaml()
    generator.create()

def run_generator_remove(args):
    generator = Generator.from_inputyaml()
    generator.remove()

def run_manager_run(args):
    assert args.run is not None, 'args.run should have a value: interactive or background.'

    match args.run:
        case 'interactive':
            os.system('./run.sh')
        case 'background':
            os.system('./rund.sh')
        case _:
            NotImplementedError('fpflow manager --run=<value> should be interactive or background.')

def run_manager_plot(args):
    assert args.plot is not None, 'args.plot should have a value: no-gui or gui.'

    match args.plot:
        case 'no-gui':
            manager = Manager().from_inputyaml()
            manager.plot(show=False)
        case 'gui':
            manager = Manager().from_inputyaml()
            manager.plot(show=True)
        case _:
            NotADirectoryError('fpflow manager --plot=<value> should be no-gui or gui.')

def fpflow():
    # Main parser. 
    parser = ArgumentParser(description='fpflow: templates, generator, manager.')
    
    # Subcommand. 
    subcommand = parser.add_subparsers(dest='subcommand', help='input/generator/manager.')

    # Input. 
    input_parser = subcommand.add_parser('input', help='')
    input_parser.add_argument('--list', action='store_true', default=None, help='List available templates.')
    input_parser.add_argument('--template', type=str, nargs='?', const='master.yaml', default=None, help='Provide template. Default is master.yaml')
    
    # Generator. 
    generator_parser = subcommand.add_parser('generator', help='')
    generator_parser.add_argument('--create', action='store_true', default=None, help='Create files.')
    generator_parser.add_argument('--remove', action='store_true', default=None, help='Remove files.')

    # Manager. 
    manager_parser = subcommand.add_parser('manager', help='')
    manager_parser.add_argument('--run', type=str, nargs='?', const='interactive', default=None, help='Can be interactive/background.')
    manager_parser.add_argument('--plot', type=str, nargs='?', const='no-gui', default=None, help='Can be no-gui/gui.')

    # Run the parser. 
    args = parser.parse_args()


    # Take actions. 
    match args.subcommand:
        case 'input':
            if args.list:
                run_input_list(args)
                return
            
            if args.template is not None:
                run_input_template(args)
                return
    
        case 'generator':
            if args.create:
                run_generator_create(args)
                return
            
            if args.remove:
                run_generator_remove(args)
                return
            
        case 'manager':
            if args.run is not None:
                run_manager_run(args)
                return 
            
            if args.plot is not None:
                run_manager_plot(args)
                return 
        
        case _:
            NotImplementedError

#endregion

#region classes
#endregion