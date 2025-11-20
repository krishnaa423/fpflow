#region modules
import numpy as np
import jmespath
from fpflow.inputs.inputyaml import InputYaml
import h5py 
import glob 
from ase.io.cube import read_cube_data
from ase.atoms import Atoms
from fpflow.structure.struct import Struct
from fpflow.structure.kpts import Kpts
from ase.units import Hartree

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BseqPpAnalysis:
    def __init__(self):
        self.inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
        struct: Struct = Struct.from_inputdict(self.inputdict)
        self.atoms: Atoms = struct.atoms[struct.struct_idx]

        self.datasets: dict = {}

    def read_evecs(self):
        files: list[str] = sorted(glob.glob('./q_*/eigenvectors.h5'))
        nQ: int = len(files)
        nS: int = jmespath.search('bse.absorption.nxct', self.inputdict)
        nk: int = len(files)
        nc: int = jmespath.search('bse.absorption.cond_bands', self.inputdict)-1
        nv: int = jmespath.search('bse.absorption.val_bands', self.inputdict)

        self.xct_evecs: np.ndarray = np.zeros([nQ, nS, nk, nc, nv], dtype='c16')
        for iQ, file in enumerate(files):
            with h5py.File(file, 'r') as f:
                xct_evec_q = f['/exciton_data/eigenvectors'][0, :, :, :, :, 0, 0] + 1j * f['/exciton_data/eigenvectors'][0, :, :, :, :, 0, 1]
                xct_evec_q /= np.linalg.norm(xct_evec_q.flatten())  # Normalize the exciton wavefunction. 
                self.xct_evecs[iQ, :, :, :, :] = xct_evec_q

        self.datasets['xct_evecs'] = self.xct_evecs

    def read_eigs(self):
        files: list[str] = sorted(glob.glob('./q_*/eigenvectors.h5'))
        nQ: int = len(files)
        nS: int = jmespath.search('bse.absorption.nxct', self.inputdict)

        self.xct_eigs: np.ndarray = np.zeros([nQ, nS], dtype='c16')
        for iQ, file in enumerate(files):
            with h5py.File(file, 'r') as f:
                self.xct_eigs[iQ, :] = f['/exciton_data/eigenvalues'][:]

        self.datasets['xct_eigs'] = self.xct_eigs / Hartree # Convert to Hartree a.u. 

    # def calc_ne(self):        # Need wannier functions, so for now, just skipping this for later. 
    #     pass

    # def calc_nh(self):
    #     pass

    def read_Qpts(self):
        Qpts: np.ndarray = np.array(Kpts.from_kgrid(
            kgrid=jmespath.search('bseq.qgrid', self.inputdict),
            is_reduced=False,
        ).bseq_qpts).astype('f8')

        self.datasets['Qpts'] = Qpts

    def read_all(self):
        self.read_evecs()
        self.read_eigs()
        # self.calc_ne()
        # self.calc_nh()
        self.read_Qpts()

    def write_all(self):
        with h5py.File('./xct.h5', 'w') as h5file:
            for key, data in self.datasets.items():
                h5file.create_dataset(key, data=data)

#endregion
