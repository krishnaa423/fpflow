#region modules
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt 
import numpy as np 
import yaml 
import h5py 
from fpflow.structure.kpath import Kpath
from fpflow.inputs.inputyaml import InputYaml
import jmespath
import os 
from fpflow.plots.plot import PlotBase, PlotType
import pandas as pd
from ase.units import Hartree, eV
import glob 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwConvergenceGwPlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.get_data_and_figures()

    def get_from_single_folder(self, dest_dir: str, color, src_dir: str=os.getcwd()):
        # Change to dest dir.
        os.chdir(dest_dir)

        # Get fermi energy. 
        tree = ET.parse('./dftelbands.xml')
        root = tree.getroot()
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree

        data = np.loadtxt('./bandstructure_inteqp.dat', skiprows=2)
        num_bands = np.unique(data[:, 1]).size
        eqp = data[:, 6].reshape(num_bands, -1).T

        self.num_bands = num_bands
        self.eqp = eqp - fermi_energy

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)
        gw_ecut: str = 'ecut' + str(jmespath.search('gw.sigma.ecut', inputdict))
        gw_bands: str = 'b' + str(jmespath.search('gw.sigma.conv_cond_bands', inputdict))
        dset_name: str = f'dset_{gw_ecut}_{gw_bands}'

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(num_bands)]

        eqp_df = pd.DataFrame(
            np.hstack([self.axis, self.eqp]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": dset_name,
                "data": eqp_df  # Store df as a single object in one row
            }
        ])

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

        for ib in range(self.num_bands):
            append_fig_df: pd.DataFrame = pd.DataFrame([{
                'fig_name': 'convqegw',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} GW Bandstructure',
                'dset_name': dset_name,
                'dset_axis_cols': 'x',        
                'dset_data_cols': [self.data_colnames[ib]],
                'color': color, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': f'{gw_ecut}_{gw_bands}' if ib==0 else None,
            }])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)


        # Change back to src dir.
        os.chdir(src_dir)

    def get_data_and_figures(self):
        dirs = glob.glob('./convergence/qe/gw/*')
        dirs.sort()

        colors = plt.cm.tab20(np.linspace(0, 1, len(dirs)))

        for idir, color in zip(dirs, colors):
            self.get_from_single_folder(idir, color)

#endregion