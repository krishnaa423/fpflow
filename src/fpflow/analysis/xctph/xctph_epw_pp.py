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
class EpwXctphPpAnalysis:
    def __init__(self):
        self.inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
        struct: Struct = Struct.from_inputdict(self.inputdict)
        self.atoms: Atoms = struct.atoms[struct.struct_idx]

        self.datasets: dict = {}

    def read_umat(self):
        with open('./struct_u.mat') as f:
            lines = f.readlines()

        # Find the first line after the header that contains three integer-like values
        shape_line = None
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    nk = int(float(parts[0]))
                    nwan1 = int(float(parts[1]))
                    nwan2 = int(float(parts[2]))
                    shape_line = (nk, nwan1, nwan2)
                    shape_idx = i
                    break
                except ValueError:
                    continue
        if shape_line is None:
            raise ValueError("Could not find valid shape line in file.")

        nk, nwan1, nwan2 = shape_line
        assert nwan1 == nwan2, "Matrix must be square"
        nwan = nwan1

        idx = shape_idx + 1  # Start after shape line
        u_matrices = []
        kpts = []
        for _ in range(nk):
            # Find k-point line (should have 3 floats)
            while idx < len(lines) and len(lines[idx].split()) < 3:
                idx += 1
            kpt = [float(x.replace('+','')) for x in lines[idx].split()[:3]]
            kpts.append(kpt)
            idx += 1
            block = []
            for _ in range(nwan * nwan):
                real_imag = lines[idx].split()
                real, imag = float(real_imag[0]), float(real_imag[1])
                block.append(complex(real, imag))
                idx += 1
            mat = np.array(block).reshape((nwan, nwan), order='C')
            u_matrices.append(mat)
            # Skip blank line(s)
            while idx < len(lines) and lines[idx].strip() == '':
                idx += 1

        self.u_matrices = np.array(u_matrices)  # shape (nk, nwan, nwan)
        self.u_kpts = np.array(kpts)  # shape (nk, 3)
        self.datasets['u_matrices'] = self.u_matrices
        self.datasets['u_kpts'] = self.u_kpts

    def read_wanwfn(self):
        real_files: list[str] = sorted(glob.glob('./struct_*r.cube'))
        imag_files: list[str] = sorted(glob.glob('./struct_*i.cube'))

        wanwfn: np.ndarray = None
        for bnd_idx, (real_file, imag_file) in enumerate(zip(real_files, imag_files)):
            # Read cube data. 
            real_data, _ = read_cube_data(real_file)
            imag_data, _ = read_cube_data(imag_file)
            complex_data = real_data + 1j * imag_data
            complex_data /= np.linalg.norm(complex_data.flatten())  # Normalize wavefunction.

            # Create the wanwfn array on the first iteration. 
            if bnd_idx == 0:
                wfnwan = np.zeros([
                    len(real_files),
                    complex_data.shape[0],
                    complex_data.shape[1],
                    complex_data.shape[2],
                ], dtype='c16')

            wfnwan[bnd_idx, :, :, :] = complex_data

        self.datasets['wanwfn'] = wfnwan

    def read_pheigs(self):
        self.pheigs: np.ndarray = np.loadtxt('./G_full_epmatq_all/phband.fmt').astype('f4')
        self.datasets['pheigs'] = self.pheigs / 2.0 # Energies are in Hartree atomic units energy. 

    def read_phevecs(self):
        files: list[str] = sorted(glob.glob('./G_full_epmatq_all/ph_eigvec_*.dat'))
        nq: int = len(files)
        nmodes: int = len(self.atoms)*3

        self.phevecs: np.ndarray = None
        for qpt_idx, file in enumerate(files):
            data = np.loadtxt(file)
            phevec: np.ndarray = (data[:, 0] + 1j * data[:, 1]).reshape(nmodes, nmodes)

            if qpt_idx == 0:
                self.phevecs = np.zeros((nq, nmodes, nmodes), dtype='c16')

            self.phevecs[qpt_idx, :, :] = phevec
        self.datasets['phevecs'] = self.phevecs

    def read_xctph(self):
        files: list[str] = sorted(glob.glob('./G_full_epmatq_all/G_full_epmatq_*.dat'))
        nQ: int = len(files)
        nS: int = jmespath.search('xctph.nxct', self.inputdict)
        nmodes: int = len(self.atoms)*3

        self.xctph: np.ndarray = None
        for qpt_idx, file in enumerate(files):
            data = np.loadtxt(file)
            xctph_q: np.ndarray = (data[:, 5] + 1j * data[:, 6]).reshape(nQ, nS, nS, nmodes)
            xctph_q = np.einsum('QsSu->sSuQ', xctph_q)      # Rearrange to convenient format. 

            # Initialize the size. 
            if qpt_idx == 0:
                self.xctph = np.zeros([nS, nS, nmodes, nQ, nQ], dtype='c16')

            self.xctph[:, :, :, :, qpt_idx] = xctph_q / Hartree   # Convert to Hartree a.u. 

        self.datasets['xctph'] = self.xctph

    def read_xctph_el(self):
        files: list[str] = sorted(glob.glob('./G_full_epmatq_el/G_full_epmatq_*.dat'))
        nQ: int = len(files)
        nS: int = jmespath.search('xctph.nxct', self.inputdict)
        nmodes: int = len(self.atoms)*3

        self.xctph: np.ndarray = None
        for qpt_idx, file in enumerate(files):
            data = np.loadtxt(file)
            xctph_q: np.ndarray = (data[:, 5] + 1j * data[:, 6]).reshape(nQ, nS, nS, nmodes)
            xctph_q = np.einsum('QsSu->sSuQ', xctph_q)      # Rearrange to convenient format. 

            # Initialize the size. 
            if qpt_idx == 0:
                self.xctph = np.zeros([nS, nS, nmodes, nQ, nQ], dtype='c16')

            self.xctph[:, :, :, :, qpt_idx] = xctph_q / Hartree   # Convert to Hartree a.u. 

        self.datasets['xctph_el'] = self.xctph

    def read_xctph_hole(self):
        files: list[str] = sorted(glob.glob('./G_full_epmatq_hole/G_full_epmatq_*.dat'))
        nQ: int = len(files)
        nS: int = jmespath.search('xctph.nxct', self.inputdict)
        nmodes: int = len(self.atoms)*3

        self.xctph: np.ndarray = None
        for qpt_idx, file in enumerate(files):
            data = np.loadtxt(file)
            xctph_q: np.ndarray = (data[:, 5] + 1j * data[:, 6]).reshape(nQ, nS, nS, nmodes)
            xctph_q = np.einsum('QsSu->sSuQ', xctph_q)      # Rearrange to convenient format. 

            # Initialize the size. 
            if qpt_idx == 0:
                self.xctph = np.zeros([nS, nS, nmodes, nQ, nQ], dtype='c16')

            self.xctph[:, :, :, :, qpt_idx] = xctph_q / Hartree   # Convert to Hartree a.u. 

        self.datasets['xctph_hole'] = self.xctph

    def read_Qpts(self):
        Qpts: np.ndarray = np.array(Kpts.from_kgrid(
            kgrid=jmespath.search('bseq.qgrid', self.inputdict),
            is_reduced=False,
        ).bseq_qpts).astype('f8')

        self.datasets['Qpts'] = Qpts

    def read_all(self):
        self.read_umat()
        self.read_wanwfn()
        self.read_pheigs()
        self.read_phevecs()
        self.read_xctph()
        self.read_xctph_el()
        self.read_xctph_hole()
        self.read_Qpts()

    def write_all(self):
        with h5py.File('./xctph.h5', 'w') as h5file:
            for key, data in self.datasets.items():
                h5file.create_dataset(key, data=data)


#endregion
