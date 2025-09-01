#region modules
import numpy as np 
import matplotlib.pyplot as plt 
import xml.etree.ElementTree as ET
import re 
import glob 
from scipy.interpolate import griddata
from typing import List, Dict, Tuple
import h5py 
import copy
from scipy import spatial
from ase.units import Ry 
import numpy as np 
from typing import List, Union
from fpflow.inputs.inputyaml import InputYaml
from fpflow.structure.kpath import Kpath
import jmespath
import os 
#endregion

#region variables
#endregion

#region functions
def unit_range(klist_in):

    tol = 1e-6
    klist = copy.copy(klist_in)

    for ik, k in enumerate(klist):
        for i, kx in enumerate(k):

            while kx < -tol:
                kx = kx + 1.0
            while kx > 1.0 -tol:
                kx = kx - 1.0

            k[i] = kx

        klist[ik,:] = k

    return klist

def find_kpt(ktargets, klist):

    klist = unit_range(klist)
    ktargets = unit_range(ktargets)

    tree = spatial.KDTree(klist)
    ik_addr = list()
    for k in ktargets:
        d, i = tree.query(k)
        if d > 1e-6:
            print('kpt not found:', k)
            i = None

        ik_addr.append(i)

    return ik_addr
#endregion

#region classes
class XctphPlot:
    def __init__(
        self,
        xctph_filename: str = 'xctph.h5',
        phbands_filename: str = 'struct.freq.gp',
        bandpathpkl_filename: str = 'bandpath.pkl',
        input_filename: str = 'fullgridflow.pkl',
        xctph_mult_factor: float=1.0e3,
        xct_Qpt_idx: int=0, # 0 based index. 
        plot_xct_state: int = 0,
    ):
        '''
        Inputs:
          xct_state: int
            A zero based index for the exciton state. Default value is 0, which indicates the lowest exciton state. 
        '''
        self.xctph_filename: str = xctph_filename
        self.phbands_filename: str = phbands_filename
        self.bandpathpkl_filename: str = bandpathpkl_filename
        self.input_filename: str = input_filename
        self.xctph_mult_factor: float = xctph_mult_factor
        self.xct_Qpt_idx: int = xct_Qpt_idx
        self.plot_xct_state: int = plot_xct_state
        self.kpath: Kpath = Kpath.from_yamlfile() 

        # Additional data created. 
        self.num_bands: int = None 
        self.phbands: np.ndarray = None
        self.inputdict: dict = InputYaml.from_yaml_file().inputdict
        self.xctph_interpolated: np.ndarray = None 

    def get_phbands_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[1]

    def get_xctph_gridded(self, xctph_values, kpts_flat):
        kgrid_size = np.array(self.inputdict['wfn']['kgrid'], dtype='i4')
        xctph_gridded = np.zeros(shape=kgrid_size)
        for x_idx in range(kgrid_size[0]):
            for y_idx in range(kgrid_size[1]):
                for z_idx in range(kgrid_size[2]):
                    kpt = np.array([
                        x_idx/kgrid_size[0],
                        y_idx/kgrid_size[1],
                        z_idx/kgrid_size[2],
                    ]).reshape(1, 3)
                    kpt_idx = find_kpt(kpt, kpts_flat)[0]
                    if kpt_idx is None or kpt_idx<0:
                        raise Exception(f'kpt_idx is not valid: {kpt_idx}')
                    xctph_gridded[x_idx, y_idx, z_idx] = xctph_values[kpt_idx]

        xctph_gridded = np.pad(xctph_gridded, pad_width=1, mode='wrap')[1:, 1:, 1:]
        kpts_gridded = np.zeros(shape=((kgrid_size[0]+1)*(kgrid_size[1]+1)*(kgrid_size[2]+1), 3))
        x, y, z = np.meshgrid(
            np.linspace(0, 1, kgrid_size[0]+1),
            np.linspace(0, 1, kgrid_size[1]+1),
            np.linspace(0, 1, kgrid_size[2]+1),
        )
        kpts_gridded[:, 0], kpts_gridded[:, 1], kpts_gridded[:, 2] = x.flatten(), y.flatten(), z.flatten()

        return kpts_gridded.reshape(-1, 3), xctph_gridded.reshape(-1, 1)

    def get_xctph_data(self):
        xctph: np.ndarray = None 
        qpts: np.ndarray = None 
        kpts = np.array(self.kpath.kpts)
        num_kpath_pts = kpts.shape[0]
        with h5py.File(self.xctph_filename, 'r') as r:
            num_xct_states = r['xctph_eh'].shape[0]
            num_modes = r['xctph_eh'].shape[3]
        self.xctph_interpolated = np.zeros(shape=(num_kpath_pts, num_modes, num_xct_states))

        for xct_state_index in range(num_xct_states):
            with h5py.File(self.xctph_filename, 'r') as r:
                xctph = np.abs(r['xctph_eh'][
                    xct_state_index,
                    xct_state_index,
                    self.xct_Qpt_idx,
                    :,
                    :
                ]*Ry)
                qpts = r['qpts'][:]
        
            kpath_pts = np.array(self.kpath.kpts)

            num_kpath_pts = kpath_pts.shape[0]
            num_modes = xctph.shape[0]
            for mode in range(num_modes):
                kpts_gridded, xctph_gridded = self.get_xctph_gridded(xctph[mode, :], qpts)
                # self.xctph_interpolated[:, mode] = griddata(qpts, xctph[mode, :], kpath_pts, method='linear')*self.xctph_mult_factor
                self.xctph_interpolated[:, mode, xct_state_index] = griddata(kpts_gridded, xctph_gridded, kpath_pts, method='linear').reshape(-1)*self.xctph_mult_factor

    def save_data(self):
        # Get some data. 
        self.get_phbands_data()
        self.get_xctph_data()
        kpts = np.array(self.kpath.kpts)

        # Save data. 
        kpt_axis = np.tile(np.arange(kpts.shape[0]).reshape(-1, 1), reps=(self.phbands.shape[1], 1))
        # Write: kpt, pheigs, xctph_interpolated. 
        num_kpts = self.phbands.shape[0]
        num_bands = self.phbands.shape[1]
        num_xct_states = self.xctph_interpolated.shape[-1]
        phbands = self.phbands.T.flatten()
        xctph = np.transpose(self.xctph_interpolated, axes=(1, 0, 2)).reshape(-1, self.xctph_interpolated.shape[-1])
        with open('plot_xctph.csv', 'w') as w:
            for row in range(kpt_axis.shape[0]):
                if row!=0 and row%num_kpts==0:
                    w.write('\n')
                w.write(f'{kpt_axis[row, 0]}, ')
                for xct_state in range(num_xct_states):
                    w.write(f'{xctph[row, xct_state]}, ')
                w.write(f'{phbands[row]} ')
                w.write('\n')

    def save_plot(self, save_filename='./plots/xctph.png', show=False, ylim=None):
        # Get some data. 
        self.get_phbands_data()
        self.get_xctph_data()
        kpts = self.kpath.kpts
        path_special_points = jmespath.search('kpath.special_points', self.inputdict)
        path_segment_npoints = jmespath.search('kpath.npoints_segment', self.inputdict)

        # Save data. 
        kpt_axis = np.tile(np.arange(kpts.shape[0]).reshape(-1, 1), reps=(self.phbands.shape[1], 1))
        # Write: kpt, pheigs, xctph_interpolated. 
        num_kpts = self.phbands.shape[0]
        num_bands = self.phbands.shape[1]
        num_xct_states = self.xctph_interpolated.shape[-1]
        phbands = self.phbands.T.flatten()
        xctph = np.transpose(self.xctph_interpolated, axes=(1, 0, 2)).reshape(-1, self.xctph_interpolated.shape[-1])
        # with open('plot_xctph.csv', 'w') as w:
        #     for row in range(kpt_axis.shape[0]):
        #         if row!=0 and row%num_kpts==0:
        #             w.write('\n')
        #         w.write(f'{kpt_axis[row, 0]}, ')
        #         for xct_state in range(num_xct_states):
        #             w.write(f'{xctph[row, xct_state]}, ')
        #         w.write(f'{phbands[row]} ')
        #         w.write('\n')

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments or total npoints. 
        ax.plot(self.phbands, color='blue')
        xaxis = np.arange(self.phbands.shape[0]).reshape(-1, 1)
        num_modes = self.phbands.shape[1]
        xaxis = np.repeat(xaxis, num_modes, axis=1)
        ax.scatter(xaxis, self.phbands, s=self.xctph_interpolated[:, :, self.plot_xct_state], color='red')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )

        # Set some labels. 
        ax.set_title(f'Phonon bands and xctph coupling for xct={self.plot_xct_state} and Qpt={self.xct_Qpt_idx}')
        ax.set_ylabel('Freq (cm-1)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        os.system('mkdir -p plots')
        fig.savefig(save_filename)
        if show: plt.show()

#endregion