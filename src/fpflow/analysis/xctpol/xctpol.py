#region modules
import numpy as np
from mpi4py import MPI
import h5py
from fpflow.analysis.dist.shape import Shape
from fpflow.analysis.logging.func_time import func_timing
from fpflow.analysis.logging.print_root import print_root
from petsc4py import PETSc
from slepc4py import SLEPc
import jmespath
from fpflow.inputs.inputyaml import InputYaml
from fpflow.structure.struct import Struct
from ase.atoms import Atoms

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class Xctpol:
    def __init__(self):
        self.inputdict: InputYaml = InputYaml.from_yaml_file('../input.yaml').inputdict
        self.temps: np.ndarray = np.array(
            jmespath.search('xctpol.temps', self.inputdict)
        )
        self.delta: float = 1e-3 
        struct: Struct = Struct.from_inputdict(self.inputdict)
        self.atoms: Atoms = struct.atoms[struct.struct_idx]

        self.use_all_xctph: bool = True # This can change for initial guess.  

        self.nQ: int = int(np.prod(jmespath.search('bseq.qgrid', self.inputdict)))
        self.nS: int = jmespath.search('xctpol.nxct', self.inputdict)
        self.nmodes: int = len(self.atoms)*3

        # Can change during loop. 
        self.current_Qpt_idx: int = 0

    def read_xctph_all(self) -> np.ndarray:
        with h5py.File('../zd_epw/xctph.h5', 'r') as xctph_file:
            xctph_all: np.ndarray = xctph_file['/xctph'][:]
        return xctph_all

    def read_xctph_hole(self) -> np.ndarray:
        with h5py.File('../zd_epw/xctph.h5', 'r') as xctph_file:
            xctph_hole: np.ndarray = xctph_file['/xctph_hole'][:]
        return xctph_hole

    def read_pheigs(self) -> np.ndarray:
        with h5py.File('../zd_epw/xctph.h5', 'r') as xctph_file:
            pheigs: np.ndarray = xctph_file['/pheigs'][:]
        return pheigs
    
    def read_xcteigs(self) -> np.ndarray:
        with h5py.File('../bseq/xct.h5', 'r') as xctph_file:
            xcteigs: np.ndarray = xctph_file['/xct_eigs'][:]
        return xcteigs
    
    def read_arrays(self):
        self.xctph_all: np.ndarray = self.read_xctph_all()
        self.xctph_hole: np.ndarray = self.read_xctph_hole()
        self.pheigs: np.ndarray = self.read_pheigs()
        self.xcteigs: np.ndarray = self.read_xcteigs()
        self.ph_occ: np.ndarray = np.zeros([self.nQ, self.nmodes], dtype='f8')

        # Set the inverse phonon eigenvalues. 
        self.inv_pheigs: np.ndarray = np.zeros([self.nQ, self.nmodes], dtype='f8')
        for iQ in range(self.nQ):
            for imode in range(self.nmodes):
                if abs(self.pheigs[iQ, imode].real) > 1e-8:
                    self.inv_pheigs[iQ, imode] = 1.0 / self.pheigs[iQ, imode].real
                else:
                    self.inv_pheigs[iQ, imode] = 0.0
        
    def setup_eps(self):
        self.vec: PETSc.Vec = PETSc.Vec().create()
        self.vec.setSizes(self.nS)
        self.vec.setType(PETSc.Vec.Type.SEQ)
        self.vec.setUp()

        self.mat: PETSc.Mat = PETSc.Mat().create()
        self.mat.setType(PETSc.Mat.Type.DENSE)
        self.mat.setSizes([self.nS, self.nS])
        self.mat.setUp()

        self.eps: SLEPc.EPS = SLEPc.EPS().create()
        self.eps.setProblemType(SLEPc.EPS.ProblemType.NHEP)
        self.eps.setOperators(self.mat)
        self.eps.setDimensions(2)
        self.eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
        self.eps.setUp()

    def init_arrays(self):
        self.tp: np.ndarray = np.zeros([self.nQ, self.nS, self.nS], dtype='c16')
        self.fm: np.ndarray = np.zeros([self.nQ, self.nS, self.nS], dtype='c16')
        self.diag: np.ndarray = np.zeros([self.nQ, self.nS, self.nS], dtype='c16')
        self.ham: np.ndarray = np.zeros([self.nQ, self.nS, self.nS], dtype='c16')
        self.eigs: np.ndarray = np.zeros([self.nQ, self.nS], dtype='c16')
        self.evecs: np.ndarray = np.zeros([self.nQ, self.nS, self.nS], dtype='c16')

        self.xctpol_eigs: np.ndarray = np.zeros([self.nQ], dtype='c16')
        self.xctpol_evecs: np.ndarray = np.ones([self.nQ, self.nS], dtype='c16') / np.sqrt(self.nS)

        self.setup_eps()

        self.prev_eig_real: float = None 

        # Set the diagonal matrix.
        for Qpt_idx in range(self.nQ):
            self.diag[Qpt_idx, :, :] = np.diag(self.xcteigs[Qpt_idx, :].real)

    def get_initial_guess(self):
        # Temporarily set to False to get initial guess.
        self.use_all_xctph = False

        self.calc_tp()
        self.assemble_ham()
        self.diagonalize_ham()

        # Set it back after getting initial guess. 
        self.use_all_xctph = True

    def calc_tp(self):
        '''
        Per qpoint index. 
        '''
        if self.use_all_xctph:
            self.tp[self.current_Qpt_idx, :, :] = np.einsum(
                'sSu, bauC, Cb, Ca, u -> sS',
                self.xctph_all[:, :, :, self.current_Qpt_idx, 0],
                self.xctph_all[:, :, :, :, 0].conj(),
                self.xctpol_evecs,
                self.xctpol_evecs.conj(),
                -2.0*self.inv_pheigs[0, :],
            )
        else:
            self.tp[self.current_Qpt_idx, :, :] = np.einsum(
                'sSu, bauC, Cb, Ca, u -> sS',
                self.xctph_hole[:, :, :, self.current_Qpt_idx, 0],
                self.xctph_hole[:, :, :, :, 0].conj(),
                self.xctpol_evecs,
                self.xctpol_evecs.conj(),
                -2.0*self.inv_pheigs[0, :],
            )

    def calc_fm(self):
        '''
        For all Qpoints. 
        '''
        pass

    def assemble_ham(self):
        self.ham[self.current_Qpt_idx, :, :] = self.tp[self.current_Qpt_idx, :, :] + self.diag[self.current_Qpt_idx, :, :]

    def diagonalize_ham(self):
        # Set the matrix for this Q-point.
        self.mat[:, :] = self.ham[self.current_Qpt_idx, :, :]
        self.mat.assemble()

        # Solve the eigenvalue problem.
        self.eps.setOperators(self.mat)
        self.eps.setUp()
        self.eps.solve()

        nconv: int = self.eps.getConverged()
        if nconv < 1:
            print_root(f'Warning: No eigenvalues converged for Q-point {self.current_Qpt_idx}.')
        
        # We only care about the lowest eigenvalue/eigenvector.
        self.xctpol_eigs[self.current_Qpt_idx] = self.eps.getEigenpair(0, self.vec)
        self.xctpol_evecs[self.current_Qpt_idx, :] = self.vec.getArray()

        # Set current eig real.
        self.current_eig_real: float = self.xctpol_eigs[self.current_Qpt_idx].real

    def step(self):
        self.calc_tp()
        self.assemble_ham()
        self.diagonalize_ham()

    def get_error(self):
        return abs(self.current_eig_real - self.prev_eig_real)

    def run_convergence(self):
        self.max_error: float = float(jmespath.search('xctpol.max_error', self.inputdict))
        self.max_steps: int = int(jmespath.search('xctpol.max_steps', self.inputdict))

        for Qpt_idx in range(self.nQ):
            self.current_Qpt_idx = Qpt_idx
            self.error: float = self.max_error + 1.0

            # Initial guess is calculated for each Qpt. 
            self.get_initial_guess()

            for step in range(self.max_steps):
                self.step()

                # Check convergence. 
                if self.prev_eig_real is None:  # Atleast have one step done. 
                    print_root(f'Qpt: {self.current_Qpt_idx}, Step {step}: Error = {self.error}')
                    self.prev_eig_real = self.current_eig_real
                    continue
                else:
                    self.error: float = self.get_error()
                    print_root(f'Qpt: {self.current_Qpt_idx}, Step {step}: Error = {self.error}')
                    if self.error < self.max_error:
                        print_root(f'Xctpol converged in {step} steps with error {self.error}.')
                        break
                    self.prev_eig_real = self.current_eig_real
    
    def run(self):
        self.read_arrays()
        self.init_arrays()
        self.run_convergence()
        self.calc_fm()

    def write(self):
        datasets: dict = {
            'xctpol_eigs': self.xctpol_eigs,
            'xctpol_evecs': self.xctpol_evecs,
            'ham': self.ham,
            'tp': self.tp,
            'diag': self.diag,
        }

        with h5py.File('./xctpol.h5', 'w') as f:
            for dset_name, data in datasets.items():
                f.create_dataset(dset_name, data=data)

#endregion