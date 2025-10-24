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
class BseAbsorptionPlot(PlotBase):
    def __init__(
        self,
        eh_filename='./absorption_eh.dat',
        noeh_filename='absorption_noeh.dat',
        outfile_prefix='bse_absorption',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.eh_filename: str = eh_filename
        self.noeh_filename: str = noeh_filename
        self.outfile_prefix: str = outfile_prefix
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        abs_eh_data = np.loadtxt(self.eh_filename, dtype='f8', skiprows=4)
        abs_noeh_data = np.loadtxt(self.noeh_filename, dtype='f8', skiprows=4)
        eh_axis = abs_eh_data[:, 0]
        eh_data = abs_eh_data[:, 1]
        noeh_axis = abs_noeh_data[:, 0]
        noeh_data = abs_noeh_data[:, 1]

        # Get name.
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_bse_absorption_eh"],
            "data": [pd.DataFrame({
                "x_eh": eh_axis,
                "y_eh": eh_data,
            })]
        })
        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_bse_absorption_noeh"],
            "data": [pd.DataFrame({
                "x_noeh": noeh_axis,
                "y_noeh": noeh_data,
            })]
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': self.outfile_prefix,
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'Energy (eV)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': r'$\epsilon_2(\omega) (arb.)$', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} BSE Absorption',
                'dset_name': 'dset_bse_absorption_noeh',
                'dset_axis_cols': 'x_noeh',        
                'dset_data_cols': ['y_noeh'],
                'color': 'blue', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': 'noeh',
            },
            {
                'fig_name': self.outfile_prefix,
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'Energy (eV)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': r'$\epsilon_2(\omega) (arb.)$', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} BSE Absorption',
                'dset_name': 'dset_bse_absorption_eh',
                'dset_axis_cols': 'x_eh',        
                'dset_data_cols': ['y_eh'],
                'color': 'red', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': 'eh',
            },
        ])

        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)            

#endregion