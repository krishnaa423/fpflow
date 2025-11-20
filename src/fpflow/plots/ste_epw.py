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
        # xctpol_100: np.ndarray = (xct_grid[:, 0] - np.array([0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0.001, 0.0015])).reshape(-1, 1)
        # xctpol_200: np.ndarray = (xct_grid[:, 0] - np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])).reshape(-1, 1)
        xctpol_300: np.ndarray = (xct_grid[:, 0] - np.array([0.401, 0.400, 0.372, 0.352, 0.311, 0.302, 0.293, 0.206])).reshape(-1, 1)
        # ste_100: np.ndarray = np.full_like(xctpol_100, 9.5)
        # ste_200: np.ndarray = np.full_like(xctpol_100, 9.4)
        ste_300: np.ndarray = np.full_like(xctpol_300, xct_grid[0, 0] - 0.487)
        self.num_xctbands: int = xct_grid.shape[1]
        self.xct_bands = interpolate.RBFInterpolator(self.Qpts, xct_grid, neighbors=8)(self.kpath.kpts)
        # self.xctpol_100 = interpolate.RBFInterpolator(self.Qpts, xctpol_100, neighbors=8)(self.kpath.kpts)
        # self.xctpol_200 = interpolate.RBFInterpolator(self.Qpts, xctpol_200, neighbors=8)(self.kpath.kpts)
        self.xctpol_300 = interpolate.RBFInterpolator(self.Qpts, xctpol_300, neighbors=8)(self.kpath.kpts)
        # self.ste_100 = interpolate.RBFInterpolator(self.Qpts, ste_100, neighbors=8)(self.kpath.kpts)
        # self.ste_200 = interpolate.RBFInterpolator(self.Qpts, ste_200, neighbors=8)(self.kpath.kpts)
        self.ste_300 = interpolate.RBFInterpolator(self.Qpts, ste_300, neighbors=8)(self.kpath.kpts)

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
        self.xct_bands_size_colnames = [f"s{i+1}" for i in range(self.num_xctbands)] 
        self.xct_bands_color_colnames = [f"c{i+1}" for i in range(self.num_xctbands)] 
        df_xct: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.xct_bands, self.AQS_xctbands*20e1, self.AQS_xctbands]),
            columns=["x"] + self.xct_bands_colnames + self.xct_bands_size_colnames + self.xct_bands_color_colnames
        )
        append_xct_bands_df = pd.DataFrame({
            "name": ["dset_xct"],
            'data': [df_xct],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_xct_bands_df], ignore_index=True)

        # Add xctpol and ste.
        df_ste: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.xctpol_300, self.ste_300]),
            columns=["x", "xctpol_300", "ste_300"],
        )
        append_ste_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_ste"],
            'data': [df_ste],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_ste_df], ignore_index=True)

        # Add ph_bands and Bqu.
        self.phbands_colnames = [f"y{i+1}" for i in range(self.num_phbands)]
        self.phbands_size_colnames = [f"s{i+1}" for i in range(self.num_phbands)]
        self.phbands_color_colnames = [f"c{i+1}" for i in range(self.num_phbands)]
        df_ph: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.phbands, self.Bqu_phbands*4, self.Bqu_phbands]),
            columns=["x"] + self.phbands_colnames + self.phbands_size_colnames + self.phbands_color_colnames
        )
        append_ph_df = pd.DataFrame({
            "name": ["dset_ph"],
            'data': [df_ph],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_ph_df], ignore_index=True)

    def set_figures(self):
        # Add xct_bands and xctph_xctbands. 
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'ste_AQS',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$A_{\mathbf{Q}S}$', 'zlim': (self.AQS_xctbands.flatten().min(), self.AQS_xctbands.flatten().max()), 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} STE Projection on Exciton Bandstructure',
                'dset_name': 'dset_xct',
                'dset_axis_cols': 'x',        
                'dset_data_cols': self.xct_bands_colnames + self.xct_bands_size_colnames + self.xct_bands_color_colnames,
                'color': None,
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
        ])
        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)
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
                    'fig_name': 'ste_xctpol',
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
            ])
            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

        # Add ste_xctpol.
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'ste_xctpol',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': (9, 12.5), 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Exciton Polaron and STE energies',
                'dset_name': 'dset_ste',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['xctpol_300'],
                'color': 'orange', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': r'Xctpol $300K$ ',
            },
            {
                'fig_name': 'ste_xctpol',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.DASHEDLINE, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': (9, 12.5), 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Exciton Polaron and STE energies',
                'dset_name': 'dset_ste',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['ste_300'],
                'color': 'orange', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': r'DiagMC STE $300K$',
            },
        ])
        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

        # Add xctph_phbands.
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'ste_Bqu',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$B_{\mathbf{q}\lambda}$', 'zlim': (self.Bqu_phbands.flatten().min(), self.Bqu_phbands.flatten().max()), 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} STE Displacement on Phonon Bandstructure',
                'dset_name': 'dset_ph',
                'dset_axis_cols': 'x',        
                'dset_data_cols': self.phbands_colnames + self.phbands_size_colnames + self.phbands_color_colnames,
                'color': None, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
        ])
        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)
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
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

#endregion