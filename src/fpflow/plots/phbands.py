#region modules
import matplotlib.pyplot as plt 
import numpy as np 
import yaml 
import h5py 
from fpflow.structure.kpath import Kpath
from fpflow.inputs.inputyaml import InputYaml
import jmespath
import os 
from fpflow.plots.common import set_common
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class PhbandsPlot:
    def __init__(
        self,
        phbands_filename: str = './struct.freq.gp',
    ):
        self.phbands_filename = phbands_filename

        self.num_bands: int = None 
        self.phbands: np.ndarray = None 
        self.kpath: Kpath = Kpath.from_yamlfile() 
        self.inputdict: dict = InputYaml.from_yaml_file().inputdict
        self.outdata_filename: str = './plots/plot_phbands.h5'
        os.system('mkdir -p ./plots')

    def get_data(self):
        data = np.loadtxt(self.phbands_filename)
        self.phbands = data[:, 1:]

        self.num_bands = self.phbands.shape[1]
        
    def save_data(self):
        # Get some data. 
        self.get_data()
        # kpts = self.kpath.get_kpts()

        with h5py.File(self.outdata_filename, 'w') as f:
            # f.create_dataset('kpts', data=kpts)
            f.create_dataset('pheigs', data=self.phbands)

    def save_plot(self, save_filename='./plots/phbands.png', show=False, ylim=None):
        # Get some data. 
        self.get_data()
        # TODO: debug kpts assignment. 
        # kpts = self.kpath.get_kpts()
        path_special_points = jmespath.search('kpath.special_points', self.inputdict)
        path_segment_npoints = jmespath.search('kpath.npoints_segment', self.inputdict)

        with h5py.File(self.outdata_filename, 'w') as f:
            # f.create_dataset('kpts', data=kpts)
            f.create_dataset('pheigs', data=self.phbands)

        plt.style.use('bmh')
        fig = plt.figure()
        ax = fig.add_subplot()

        # Set xaxis based on segments.
        ax.plot(self.phbands, color='blue')
        ax.yaxis.grid(False)  
        ax.set_xticks(
            ticks=np.arange(len(path_special_points))*path_segment_npoints,
            labels=path_special_points,
        )

        # Set some labels. 
        ax.set_title('Phonon Bandstructure')
        ax.set_ylabel('Freq (cm-1)')
        if ylim: ax.set_ylim(bottom=ylim[0], top=ylim[1])
        os.system('mkdir -p plots')
        fig.savefig(save_filename)
        if show: plt.show()

class PhonopyPlot(PhbandsPlot):
    def __init__(self, **kwargs):
        super().__init__(phbands_filename='band.yaml', **kwargs)
        self.outdata_filename: str = './plots/plot_phonopy.h5'

    def get_data(self):
        with open(self.phbands_filename) as f: data = yaml.safe_load(f)

        nk = len(data['phonon'])
        nb = len(data['phonon'][0]['band'])

        # fill phbands
        self.phbands = np.zeros(shape=(nk, nb), dtype='f8')
        for (k, b), value in np.ndenumerate(self.phbands):
            self.phbands[k, b] = data['phonon'][k]['band'][b]['frequency']*33.356        # Factor in cm^{-1}
 
        self.num_bands = self.phbands.shape[1]

#endregion