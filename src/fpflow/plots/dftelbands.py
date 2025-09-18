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

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class DftelbandsPlot(PlotBase):
    def __init__(
        self,
        infile='./dftelbands.xml',
        outfile_prefix='dftelbands',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.infilename: str = infile
        self.outfile_prefix: str = outfile_prefix
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        tree = ET.parse(self.infilename)
        root = tree.getroot()

        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
        num_kpts = len(eig_nodes)
        num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - fermi_energy

        self.dft_eigs = dft_eigs
        self.numbands = dft_eigs.shape[1]

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file().inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(self.dft_eigs.shape[1])]

        # Build dataframe with x + y's
        df = pd.DataFrame(
            np.hstack([self.axis, self.dft_eigs]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_dftelbands"],
            "data": [df]  # Store df as a single object in one row
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        for ib in range(self.numbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([{
                'fig_name': self.outfile_prefix,
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': None, 'xlim': None, 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} DFT Bandstructure',
                'dset_name': 'dset_dftelbands',
                'dset_axis_cols': 'x',        
                'dset_data_cols': [self.data_colnames[ib]],
                'color': 'blue', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            }])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)            

#endregion