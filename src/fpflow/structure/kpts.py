#region modules
import numpy as np 
from typing import Iterable
import itertools 
from fpflow.structure.struct import Struct
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.io.logging import get_logger
from ase import Atoms 
from fpflow.structure.struct import Struct
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class Kpts:
    def __init__(self, *args, **kwargs):
        self.is_grid: bool = False
        self.is_reduced: bool = False
        self.kpts: Iterable = None 
        self.kgrid: Iterable = None 
        self.qshift: Iterable = None 

        for key, value in kwargs.items():
            setattr(self, key, value)

    def populate_kpts(self, **kwargs):
        self.kpts = np.array(list(itertools.product(range(self.kgrid[0]), range(self.kgrid[1]), range(self.kgrid[2])))) / np.array(self.kgrid)

    @classmethod
    def from_kgrid(cls, kgrid: Iterable, qshift: Iterable=[0.0, 0.0, 0.0], is_reduced: bool = False, **kwargs):
        class_call = None 
        
        match is_reduced:
            case False:
                class_call = Kpts
            case True:
                class_call = IbzKpts

        result = class_call(
            kgrid=kgrid,
            qshift=qshift,
            is_reduced=is_reduced,
            is_grid=True,
        ) 
        result.populate_kpts(**kwargs)

        return result 
    
    @property
    def nkpt(self):
        if isinstance(self.kpts, list):
            return len(self.kpts)
        elif isinstance(self.kpts, np.ndarray):
            return self.kpts.shape[0]
        else:
            return 0

    @property
    def wfn_kpts(self):
        kpts = self.kpts.tolist()
        for row in kpts: row.append(1.0)

        return kpts
    
    @property
    def wfnq_kpts(self):
        kpts = self.kpts.tolist()
        for row in kpts: row.append(1.0)

        for row in kpts:
            row[0] += self.qshift[0]
            row[1] += self.qshift[1]
            row[2] += self.qshift[2]

        return kpts

    @property
    def epsilon_kpts(self):
        kpts = self.kpts.tolist()
        for row in kpts: row.append(1.0); row.append(0)

        kpts[0][0] = self.qshift[0]
        kpts[0][1] = self.qshift[1]
        kpts[0][2] = self.qshift[2]
        kpts[0][4] = 1

        return kpts

    @property
    def sigma_kpts(self):
        kpts = self.kpts.tolist()
        for row in kpts: row.append(1.0)

        return kpts 
    
    @property
    def bseq_qpts(self):
        kpts = self.kpts.tolist()

        return kpts

class IbzKpts(Kpts):
    def _ir_kpoint_groups(self, kgrid, atoms, shift=(0, 0, 0), time_reversal=True, symprec=1e-5):
        """
        Return symmetry-equivalence classes of k-points on a uniform grid.

        Parameters
        ----------
        kgrid : tuple(int, int, int)
            Mesh shape (nkx, nky, nkz).
        atoms : ase.Atoms
            Structure with cell, pbc, numbers, and scaled positions.
        shift : tuple(int, int, int), optional
            Monkhorst-Pack half-grid shifts per axis (0 or 1). (0,0,0) = Γ-centered.
        time_reversal : bool, optional
            Whether to use time-reversal symmetry in the reduction.
        symprec : float, optional
            Tolerance for symmetry finding (Å).
        angle_tolerance : float, optional
            Angle tolerance (deg). Use -1.0 to let spglib choose.

        Returns
        -------
        groups : list[list[np.ndarray]]
            Each inner list contains the fractional k-points (in reduced coordinates,
            w.r.t. reciprocal lattice basis, wrapped into [0,1)) that are symmetry-equivalent.
            The first element of each list is the representative (irreducible) k-point.
        reps : np.ndarray, shape (nir, 3)
            Array of the representative (irreducible) k-points (same fractional convention).
        weights : np.ndarray, shape (nir,)
            Multiplicities (number of grid points folded into each representative).
        """
        try:
            import spglib as spg
        except Exception as e:
            raise ImportError("This function requires the 'spglib' Python package.") from e

        # Build spglib 'cell' tuple: (lattice (3x3), positions(fractional Nx3), numbers(N))
        lattice = np.array(atoms.cell.array, dtype=float)
        if not np.all(atoms.pbc):
            raise ValueError("Atoms must be periodic in all three directions for k-point symmetry reduction.")
        scaled_pos = atoms.get_scaled_positions()
        numbers = atoms.get_atomic_numbers()
        cell = (lattice, scaled_pos, numbers)

        mesh = np.asarray(kgrid, dtype=int)
        is_shift = np.asarray(shift, dtype=int)

        # Call spglib; handle both older and newer return conventions.
        try:
            # Older spglib: returns mapping (len=prod(mesh)) and grid addresses (N,3)
            mapping, grid = spg.get_ir_reciprocal_mesh(
                mesh=mesh,
                cell=cell,
                is_shift=is_shift,
                is_time_reversal=time_reversal,
                symprec=symprec,
            )
        except TypeError:
            # Newer spglib may return a dict if return_type changed; normalize to mapping+grid.
            out = spg.get_ir_reciprocal_mesh(
                mesh=mesh,
                cell=cell,
                is_shift=is_shift,
                is_time_reversal=time_reversal,
                symprec=symprec,
            )
            if isinstance(out, dict):
                mapping = out["mapping"]
                grid = out["grid_address"]
            else:
                raise

        mapping = np.asarray(mapping, dtype=int)
        grid = np.asarray(grid, dtype=int)  # integer addresses in reciprocal lattice coordinates

        # Convert integer grid addresses to fractional coordinates in reciprocal basis.
        # spglib uses units of 1/mesh step; wrap into [0,1).
        frac = (grid / mesh[None, :]) % 1.0

        # Group by representative index given in 'mapping'
        reps_idx, inverse = np.unique(mapping, return_inverse=True)
        nir = reps_idx.size
        groups = [[] for _ in range(nir)]
        for i, g in enumerate(inverse):
            groups[g].append(frac[i])

        # Convert lists to arrays and pick representatives (first in each group)
        groups = [ [np.asarray(k) for k in group] for group in groups ]
        reps = np.array([groups[i][0] for i in range(nir)])
        weights = np.array([len(g) for g in groups], dtype=int)

        return groups, reps, weights

    def populate_kpts(self, **kwargs):
        struct: Struct = Struct.from_yaml_file('./input.yaml')
        atoms: Atoms = struct.atoms[struct.struct_idx]

        # Get symmetry. 
        _, self.kpts, _ = self._ir_kpoint_groups(self.kgrid, atoms)
    
#endregion