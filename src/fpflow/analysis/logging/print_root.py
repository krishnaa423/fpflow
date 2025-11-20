#region modules
from mpi4py import MPI

#endregion

#region variables
#endregion

#region functions
def print_root(message: str, flush: bool = True) -> None:
    """Print message only from the root process (rank 0) in an MPI environment.

    Args:
        message (str): The message to print.
        flush (bool, optional): Whether to flush the output buffer. Defaults to True.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        print(message, flush=flush)

#endregion

#region classes
#endregion
