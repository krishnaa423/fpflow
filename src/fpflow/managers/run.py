#region modules
import subprocess
import time 
import os 
#endregion

#region variables
#endregion

#region functions
def subprocess_run(script, **kwargs):
    '''
    Run each script and write out some logging info. 
    '''

    current_dir = os.getcwd()
    
    # dest_dir. 
    dest_dir: str = None 
    if 'dest_dir' in kwargs.keys(): dest_dir = kwargs['dest_dir']
    if dest_dir is not None: os.chdir(dest_dir)

    #total_time.
    total_time: float = 0.0
    if 'total_time' in kwargs.keys(): total_time = kwargs['total_time']

    # Run and time. 
    start_time = time.time()
    print(f'Starting {script}.', flush=True)
    ps_result = subprocess.run(f'{script}')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    total_time += elapsed_time

    # Switch back directories. 
    os.chdir(current_dir)

    if ps_result.returncode == 0:  # Success.
        print(f'Done with {script} in {elapsed_time:15.10f} seconds in dir: {dest_dir if dest_dir is not None else current_dir}.\n\n', flush=True)
    else:               # Fail.
        print(f'Error finishing: {script}. Exited with code {ps_result.returncode}. Time elapsed is {elapsed_time:15.10f} seconds.\n\n', flush=True)
        print(f'Total time for workflow run in {total_time:15.10f} seconds.\n', flush=True)
        os._exit(ps_result.returncode)

    return total_time

#endregion

#region classes
#endregion