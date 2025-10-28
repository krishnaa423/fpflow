#region modules
from fpflow.inputs.inputyaml import InputYaml
import numpy as np
import re 
import jmespath
import h5py 
from fpflow.structure.kpath import Kpath
from scipy.interpolate import interp1d

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class WannierQeAnalysis:
    def __init__(self):
        self.inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict

    def read_eig(self, filename: str = './wan_band.dat', kpt_label_file: str = './wan_band.labelinfo.dat'):
        # Read the wannier kpath.
        self.wannier_kpt_idxs: np.ndarray = np.loadtxt(kpt_label_file, usecols=1) - 1
        self.dft_kpt_idxs: np.ndarray = np.arange(
            len(jmespath.search('kpath.special_points', self.inputdict))
        )*jmespath.search('kpath.npoints_segment', self.inputdict)

        # Read the dftelbands kpath.
        kpath: Kpath = Kpath.from_yamlfile('../input.yaml')
        self.dft_kpts = kpath.kpts
        self.dft_axis, self.xticks, self.xtick_labels = kpath.even_spaced_axis

        # Read the wannier bands data.
        arrays = []
        with open(filename) as f:
            content = f.read()
        blocks = [block for block in re.split(r'\n\s*\n', content) if block.strip()]
        for block in blocks:
            lines = [line for line in block.strip().split('\n') if line.strip()]
            arr = np.array([list(map(float, line.split())) for line in lines])
            arrays.append(arr)

        self.eigs = np.stack(arrays, axis=1)[:, :, 1]

        # Interpolate.
        interp = interp1d(self.wannier_kpt_idxs, self.dft_kpt_idxs, kind='linear', fill_value='extrapolate')
        self.wan_axis = interp(np.arange(self.eigs.shape[0]))

    def read_kgrid(self):
        self.kgrid: list = jmespath.search('wannier.kgrid', self.inputdict)

    def read_hr(self, filename: str = './wan_hr.dat'):
        with open(filename) as f:
            lines = f.readlines()

        # Parse nwan and nrpts
        nwan = int(lines[1].strip())
        nrpts = int(lines[2].strip())

        # Parse degeneracies (nrpts integers, 15 per line)
        deg_lines = []
        idx = 3
        while len(deg_lines) < nrpts:
            deg_lines += [int(x) for x in lines[idx].split()]
            idx += 1

        # Now parse the Hamiltonian matrix elements
        # Each line: R1 R2 R3 m n real imag
        hr_data = []
        rpts = []
        for i in range(nwan * nwan * nrpts):
            parts = lines[idx].split()
            R = [int(parts[0]), int(parts[1]), int(parts[2])]
            m = int(parts[3]) - 1  # 1-based to 0-based
            n = int(parts[4]) - 1
            real = float(parts[5])
            imag = float(parts[6])
            hr_data.append((i // (nwan * nwan), m, n, complex(real, imag)))
            if m == 0 and n == 0:
                rpts.append(R)
            idx += 1

        # Build (nrpts, nwan, nwan) array
        hr_matrices = np.zeros((nrpts, nwan, nwan), dtype=complex)
        for r, m, n, val in hr_data:
            hr_matrices[r, m, n] = val

        self.hr_matrices = hr_matrices
        self.hr_rpts = np.array(rpts)  # shape (nrpts, 3)

    def read_u_matrix(self, filename: str = './wan_u.mat'):
        with open(filename) as f:
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
            mat = np.array(block).reshape((nwan, nwan), order='F')
            u_matrices.append(mat)
            # Skip blank line(s)
            while idx < len(lines) and lines[idx].strip() == '':
                idx += 1

        self.u_matrices = np.array(u_matrices)  # shape (nk, nwan, nwan)
        self.u_kpts = np.array(kpts)  # shape (nk, 3)

    def read_all(self):
        self.read_eig()
        self.read_kgrid()
        self.read_hr()
        self.read_u_matrix()

    def write(self):
        datasets = {
            'eigs': self.eigs,
            'wan_axis': self.wan_axis,
            'u': self.u_matrices,
            'u_kpts': self.u_kpts,
            'hr': self.hr_matrices,
            'hr_R': self.hr_rpts, 
        }

        with h5py.File('./wannier.h5', 'w') as hf:
            for name, data in datasets.items():
                hf.create_dataset(name, data=data)

#endregion