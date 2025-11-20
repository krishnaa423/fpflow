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
class EpwXctphPlot(PlotBase):
    '''
    It will plot:
    - Exciton bandstructure along kpath.
    - Xctph projected phonon bandstructure.
    - Xctph projected exciton bandstructure.
    '''
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

        # Get phbands.
        data = np.loadtxt('./dfpt/struct.freq.gp')
        self.phbands = data[:, 1:]
        self.phbands *= 0.123984   # Convert to meV. 1 cm-1 = 0.123984 meV
        self.num_phbands: int = self.phbands.shape[1]

        # Get xctgrid. 
        with h5py.File('./bseq/xct.h5', 'r') as xct_file:
            xct_grid: np.ndarray = (xct_file['/xct_eigs'][:]).real * Hartree  # Convert to eV.
        self.num_xctbands: int = xct_grid.shape[1]

        # Interpolate xct_grid to kpath to get xct_bands.
        self.xct_bands = interpolate.RBFInterpolator(self.Qpts, xct_grid, neighbors=8)(self.kpath.kpts)

        # Get xctph_grid.
        with h5py.File('./zd_epw/xctph.h5', 'r') as xctph_file:
            self.xctph_grid: np.ndarray = xctph_file['/xctph'][:] * Hartree * 1000 # Convert to meV.

        # Interpolate xctph_grid to kpath to get xctph_phbands and xctph_xctbands.
        # self.xctph_phgrid = np.sqrt(np.einsum('sSuQq->qu', np.abs(self.xctph_grid)**2))
        # self.xctph_phgrid /= self.num_xctbands*self.num_xctbands*self.num_Qpts # Normalize.
        self.xctph_phgrid = np.einsum('uq->qu', np.abs(self.xctph_grid[0, 0, :, 0, :]))
        self.xctph_phbands = interpolate.RBFInterpolator(self.Qpts, self.xctph_phgrid, neighbors=8)(self.kpath.kpts)

        self.xctph_xctgrid = np.sqrt(np.einsum('sSuQq->QS', np.abs(self.xctph_grid)**2))
        self.xctph_xctgrid /= self.num_Qpts*self.num_phbands # Normalize.
        self.xctph_xctbands = interpolate.RBFInterpolator(self.Qpts, self.xctph_xctgrid, neighbors=8)(self.kpath.kpts)

        # Add xct_bands and xctph_xctbands. 
        self.xct_bands_colnames = [f"y{i+1}" for i in range(self.num_xctbands)]
        self.xct_bands_size_colnames = [f"s{i+1}" for i in range(self.num_xctbands)] 
        self.xct_bands_color_colnames = [f"c{i+1}" for i in range(self.num_xctbands)] 
        df_xct: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.xct_bands, self.xctph_xctbands, self.xctph_xctbands]),
            columns=["x"] + self.xct_bands_colnames + self.xct_bands_size_colnames + self.xct_bands_color_colnames
        )
        append_xct_bands_df = pd.DataFrame({
            "name": ["dset_xct"],
            'data': [df_xct],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_xct_bands_df], ignore_index=True)

        # Add xctph_phbands.
        self.phbands_colnames = [f"y{i+1}" for i in range(self.num_phbands)]
        self.phbands_size_colnames = [f"s{i+1}" for i in range(self.num_phbands)]
        self.phbands_color_colnames = [f"c{i+1}" for i in range(self.num_phbands)]
        df_ph: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.phbands, self.xctph_phbands*2e-1, self.xctph_phbands]),
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
                'fig_name': 'xctph_epw_xctph_xctbands',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$\mathcal{G}$ (meV)', 'zlim': (self.xctph_xctbands.min(), self.xctph_xctbands.max()), 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Xctph Projection on Exciton Bandstructure',
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
                    'fig_name': 'xctph_epw_xctbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Exciton Bandstructure',
                    'dset_name': 'dset_xct',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.xct_bands_colnames[band_idx]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'xctph_epw_xctph_xctbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Xctph Projection on Exciton Bandstructure',
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

        # Add xctph_phbands.
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'xctph_epw_xctph_phbands',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$\mathcal{G}$ (meV)', 'zlim': (self.xctph_phbands.min(), self.xctph_phbands.max()), 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Xctph Projection on Phonon Bandstructure',
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
                    'fig_name': 'xctph_epw_xctph_phbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Xctph Projection on Phonon Bandstructure',
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