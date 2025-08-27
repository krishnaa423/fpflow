#region modules
from argparse import ArgumentParser 
#endregion

#region variables
#endregion

#region functions
def fpflow():
    parser = ArgumentParser(description="Run fpflow tasks")
    parser.add_argument("task", type=str, help="Task to run")
    args = parser.parse_args()

    print(f"Running task: {args.task}")
#endregion

#region classes
#endregion