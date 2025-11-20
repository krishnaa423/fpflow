#region modules
from typing import Sequence
import numpy as np
from mpi4py import MPI

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class Shape:
    def __init__(self, shape: Sequence[int], world_size: int = 1, rank: int = 0):
        self.shape: np.ndarray = np.array(shape)
        self.world_size: int = world_size
        self.rank: int = rank

    @classmethod
    def from_world_comm(cls, shape: Sequence[int]) -> 'Shape':
        comm = MPI.COMM_WORLD
        return cls(shape, comm.Get_size(), comm.Get_rank())

    @property
    def size(self) -> int:
        return int(np.prod(self.shape))

    @property
    def local_size(self) -> int:
        local = self.size // self.world_size + (1 if self.rank < self.size % self.world_size else 0)
        return local

    @property
    def owner_range(self) -> tuple[int, int]:
        start: int = (self.size // self.world_size) * self.rank + min(self.rank, self.size % self.world_size)
        end: int = start + self.local_size
        return (start, end)
    

#endregion
