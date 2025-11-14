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
from fpflow.structure.kpts import Kpts
from scipy import interpolate

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwStePlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get name.
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Get Qpts.
        self.Qpts: np.ndarray = np.array(Kpts.from_kgrid(
            kgrid=[
                jmespath.search('bseq.qgrid[0]', self.inputdict),
                jmespath.search('bseq.qgrid[1]', self.inputdict),
                jmespath.search('bseq.qgrid[2]', self.inputdict),
            ],
            is_reduced=False,
        ).bseq_qpts)
        self.num_Qpts: int = self.Qpts.shape[0]

        # Get kpath.
        self.kpath: Kpath = Kpath.from_yamlfile()
        self.xaxis, self.xticks, self.xtick_labels = self.kpath.even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get phonon bands. 
        data = np.loadtxt('./dfpt/struct.freq.gp')
        self.phbands = data[:, 1:]
        self.phbands *= 0.123984   # Convert to meV. 1 cm-1 = 0.123984 meV
        self.num_phbands: int = self.phbands.shape[1]

        # get exciton bands. 
        data = np.loadtxt('./zd_epw/G_full_epmatq/exband.fmt')
        xct_grid = data[:, :]
        self.num_xctbands: int = xct_grid.shape[1]
        self.xct_bands = interpolate.RBFInterpolator(self.Qpts, xct_grid, neighbors=8)(self.kpath.kpts)

        # Get AQS on grid and interp on kpath.
        data: np.ndarray = np.loadtxt('./zd_epw/Asq.plrn', skiprows=1)
        self.AQS_grid: np.ndarray = (np.abs(data[:, 2] + 1j* data[:, 3])**2).reshape(self.num_Qpts, self.num_xctbands)
        self.AQS_xctbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.AQS_grid, neighbors=8)(self.kpath.kpts)

        # Get Bqu on grid and interp on kpath. 
        data: np.ndarray = np.loadtxt('./zd_epw/bmat.plrn', skiprows=1)
        self.Bqu_grid: np.ndarray = data[:, 5].reshape(self.num_Qpts, self.num_phbands)
        self.Bqu_phbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.Bqu_grid, neighbors=8)(self.kpath.kpts)

        # Add xct_bands and AQS. 
        self.xct_bands_colnames = [f"y{i+1}" for i in range(self.num_xctbands)]
        self.xct_bands_weight_colnames = [f"s{i+1}" for i in range(self.num_xctbands)] 
        df_xct: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.xct_bands, self.AQS_xctbands*1e2]),
            columns=["x"] + self.xct_bands_colnames + self.xct_bands_weight_colnames
        )
        append_xct_bands_df = pd.DataFrame({
            "name": ["dset_xct"],
            'data': [df_xct],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_xct_bands_df], ignore_index=True)

        # Add ph_bands and Bqu.
        self.phbands_colnames = [f"y{i+1}" for i in range(self.num_phbands)]
        self.phbands_weight_colnames = [f"s{i+1}" for i in range(self.num_phbands)]
        df_ph: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.phbands, self.Bqu_phbands*1e1]),
            columns=["x"] + self.phbands_colnames + self.phbands_weight_colnames
        )
        append_ph_df = pd.DataFrame({
            "name": ["dset_ph"],
            'data': [df_ph],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_ph_df], ignore_index=True)

    def set_figures(self):
        # Add xct_bands and xctph_xctbands. 
        for band_idx in range(self.num_xctbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': 'ste_AQS',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Projection on Exciton Bandstructure',
                    'dset_name': 'dset_xct',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.xct_bands_colnames[band_idx]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'ste_AQS',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Projection on Exciton Bandstructure',
                    'dset_name': 'dset_xct',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.xct_bands_colnames[band_idx], self.xct_bands_weight_colnames[band_idx]],
                    'color': 'red', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

        # Add xctph_phbands.
        for band_idx in range(self.num_phbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': 'ste_Bqu',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Displacement Projection on Phonon Bandstructure',
                    'dset_name': 'dset_ph',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.phbands_colnames[band_idx]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'ste_Bqu',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Displacement on Phonon Bandstructure',
                    'dset_name': 'dset_ph',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.phbands_colnames[band_idx], self.phbands_weight_colnames[band_idx]],
                    'color': 'red', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

#endregion