#region modules
import time 
from mpi4py import MPI
#endregion

#region variables
#endregion

#region functions
# Wrapper to print timing for function calls on master rank only
def func_timing(func):
    def wrapper(*args, **kwargs):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank() if MPI is not None else 0
        if rank == 0:
            print(f"[Start: {func.__qualname__}]", flush=True)
            start = time.time()
            result = func(*args, **kwargs)
            end = time.time()
            print(f"[Done: {func.__qualname__}]  Took {end - start:.6f} seconds", flush=True)
            return result
        else:
            return func(*args, **kwargs)
    return wrapper

#endregion

#region classes
#endregion
