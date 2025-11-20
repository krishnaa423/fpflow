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
class EpwPolPlot(PlotBase):
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

        # Get dftbands.
        tree = ET.parse('./dftelbands/dftelbands.xml')
        root = tree.getroot()
        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
        num_kpts = len(eig_nodes)
        num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - fermi_energy
        self.elbands = dft_eigs[:, :8]
        self.num_elbands = self.elbands.shape[1]

        # Get phbands.
        data = np.loadtxt('./dfpt/struct.freq.gp')
        self.phbands = data[:, 1:]
        self.phbands *= 0.123984   # Convert to meV. 1 cm-1 = 0.123984 meV
        self.num_phbands: int = self.phbands.shape[1]

        # Get Akn.
        data: np.ndarray = np.loadtxt('./pol_epw/elec_Ank.plrn', skiprows=1)
        self.elec_Akn_grid: np.ndarray = data[:, -1].reshape(self.num_Qpts, self.num_elbands)
        self.elec_Akn_grid: np.ndarray = self.elec_Akn_grid[:, ::-1]
        self.elec_Akn_elbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.elec_Akn_grid, neighbors=8)(self.kpath.kpts)
        data: np.ndarray = np.loadtxt('./pol_epw/hole_Ank.plrn', skiprows=1)
        self.hole_Akn_grid: np.ndarray = data[:, -1].reshape(self.num_Qpts, self.num_elbands)
        self.hole_Akn_grid: np.ndarray = self.hole_Akn_grid[:, ::-1]
        self.hole_Akn_elbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.hole_Akn_grid, neighbors=8)(self.kpath.kpts)

        # Get Bqu.
        data: np.ndarray = np.loadtxt('./pol_epw/elec_Bmat.plrn', skiprows=1)
        self.elec_Bqu_grid: np.ndarray = data[:, -1].reshape(self.num_Qpts, self.num_phbands)
        self.elec_Bqu_elbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.elec_Bqu_grid, neighbors=8)(self.kpath.kpts)
        data: np.ndarray = np.loadtxt('./pol_epw/hole_Bmat.plrn', skiprows=1)
        self.hole_Bqu_grid: np.ndarray = data[:, -1].reshape(self.num_Qpts, self.num_phbands)
        self.hole_Bqu_elbands: np.ndarray = interpolate.RBFInterpolator(self.Qpts, self.hole_Bqu_grid, neighbors=8)(self.kpath.kpts)

        # Add elbands and Akn. 
        self.elbands_colnames = [f"y{i+1}" for i in range(self.num_elbands)]
        self.elbands_elec_weight_colnames = [f"se{i+1}" for i in range(self.num_elbands)] 
        self.elbands_elec_color_colnames = [f"ce{i+1}" for i in range(self.num_elbands)] 
        self.elbands_hole_weight_colnames = [f"sh{i+1}" for i in range(self.num_elbands)] 
        self.elbands_hole_color_colnames = [f"ch{i+1}" for i in range(self.num_elbands)] 
        df_el: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.elbands, self.elec_Akn_elbands, self.hole_Akn_elbands, self.elec_Akn_elbands, self.hole_Akn_elbands]),
            columns=["x"] + self.elbands_colnames + self.elbands_elec_weight_colnames + self.elbands_hole_weight_colnames + self.elbands_elec_color_colnames + self.elbands_hole_color_colnames
        )
        append_el_bands_df = pd.DataFrame({
            "name": ["dset_el"],
            'data': [df_el],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_el_bands_df], ignore_index=True)

        # Add ph_bands and Bqu.
        self.phbands_colnames = [f"y{i+1}" for i in range(self.num_phbands)]
        self.phbands_elec_weight_colnames = [f"se{i+1}" for i in range(self.num_phbands)]
        self.phbands_elec_color_colnames = [f"ce{i+1}" for i in range(self.num_phbands)]
        self.phbands_hole_weight_colnames = [f"sh{i+1}" for i in range(self.num_phbands)]
        self.phbands_hole_color_colnames = [f"ch{i+1}" for i in range(self.num_phbands)]
        df_ph: pd.DataFrame = pd.DataFrame(
            np.hstack([self.axis, self.phbands, self.elec_Bqu_elbands, self.hole_Bqu_elbands, self.elec_Bqu_elbands, self.hole_Bqu_elbands/120*5.7]),
            columns=["x"] + self.phbands_colnames + self.phbands_elec_weight_colnames + self.phbands_hole_weight_colnames + self.phbands_elec_color_colnames + self.phbands_hole_color_colnames
        )
        append_ph_df = pd.DataFrame({
            "name": ["dset_ph"],
            'data': [df_ph],
        })
        self.dsets_df = pd.concat([self.dsets_df, append_ph_df], ignore_index=True)

    def set_figures(self):
        # Add elbands and phbands. 
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'pol_elec_Akn',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$A_{\mathbf{k}n}$', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Electron Polaron Projection on DFT Bandstructure',
                'dset_name': 'dset_el',
                'dset_axis_cols': 'x',        
                'dset_data_cols': self.elbands_colnames + self.elbands_elec_weight_colnames + self.elbands_elec_color_colnames,
                'color': None, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
            {
                'fig_name': 'pol_hole_Akn',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$A_{\mathbf{k}n}$', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Hole Polaron Projection on DFT Bandstructure',
                'dset_name': 'dset_el',
                'dset_axis_cols': 'x',        
                'dset_data_cols': self.elbands_colnames + self.elbands_hole_weight_colnames + self.elbands_hole_color_colnames,
                'color': None, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
        ])
        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)
        for band_idx in range(self.num_elbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': 'pol_elec_Akn',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Electron Polaron Projection on DFT Bandstructure',
                    'dset_name': 'dset_el',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.elbands_colnames[band_idx]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                
                {
                    'fig_name': 'pol_hole_Akn',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Hole Polaron Projection on DFT Bandstructure',
                    'dset_name': 'dset_el',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.elbands_colnames[band_idx]],
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
                'fig_name': 'pol_elec_Bqu',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': r'$B_{\mathbf{q}\lambda}$', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Electron Polaron Displacement on Phonon Bandstructure',
                'dset_name': 'dset_ph',
                'dset_axis_cols': 'x',        
                'dset_data_cols': self.phbands_colnames + self.phbands_elec_weight_colnames + self.phbands_elec_color_colnames,
                'color': None, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
            {
                    'fig_name': 'pol_hole_Bqu',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTERSIZECOLOR, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': r'$\Delta \omega_{\mathbf{q}\lambda}$ $cm^{-1}$', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Phonon frequency shift',
                    'dset_name': 'dset_ph',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': self.phbands_colnames + self.phbands_hole_weight_colnames + self.phbands_hole_color_colnames,
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
                    'fig_name': 'pol_elec_Bqu',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Electron Polaron Displacement on Phonon Bandstructure',
                    'dset_name': 'dset_ph',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.phbands_colnames[band_idx]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                
                {
                    'fig_name': 'pol_hole_Bqu',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Phonon frequency shift',
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