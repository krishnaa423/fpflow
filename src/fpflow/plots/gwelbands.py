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
class GwelbandsPlot(PlotBase):
    def __init__(
        self,
        infile='./bandstructure_inteqp.dat',
        outfile_prefix='gwelbands',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.infilename: str = infile
        self.outfile_prefix: str = outfile_prefix
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get fermi energy. 
        tree = ET.parse('./dftelbands.xml')
        root = tree.getroot()
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree

        data = np.loadtxt(self.infilename, skiprows=2)
        num_bands = np.unique(data[:, 1]).size
        emf = data[:, 5].reshape(num_bands, -1).T
        eqp = data[:, 6].reshape(num_bands, -1).T

        self.num_bands = num_bands
        self.emf = emf - fermi_energy
        self.eqp = eqp - fermi_energy

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file().inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(num_bands)]

        # Build dataframe with x + y's
        emf_dft = pd.DataFrame(
            np.hstack([self.axis, self.emf]),
            columns=["x"] + self.data_colnames
        )
        eqp_df = pd.DataFrame(
            np.hstack([self.axis, self.eqp]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": "dset_dftelbands",
                "data": emf_dft  # Store df as a single object in one row
            },
            {
                "name": "dset_gwelbands",
                "data": eqp_df  # Store df as a single object in one row
            }
        ])

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        for ib in range(self.num_bands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': self.outfile_prefix,
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': None, 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} DFT+GW Bandstructure',
                    'dset_name': 'dset_dftelbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.data_colnames[ib]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': 'DFT' if ib==0 else None,
                },
                {
                    'fig_name': self.outfile_prefix,
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': None, 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} DFT+GW Bandstructure',
                    'dset_name': 'dset_gwelbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.data_colnames[ib]],
                    'color': 'green', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': 'GW' if ib==0 else None,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)            

#endregion