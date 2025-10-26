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
from fpflow.analysis.wannierqe import WannierQeAnalysis

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeWannierPlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        
        self.get_data()
        self.set_figures()

    def add_dft_data(self):
        # Read dft eigs.
        tree = ET.parse('./dftelbands/dftelbands.xml')
        root = tree.getroot()

        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        self.fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
        num_kpts = len(eig_nodes)
        self.dft_num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, self.dft_num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - self.fermi_energy

        self.dft_eigs = dft_eigs
        self.num_bands = dft_eigs.shape[1]

        # Create column names: y1, y2, ..., yN. 
        self.dft_data_colnames = [f"y{i+1}" for i in range(self.dft_num_bands)]

        # Build dataframe with x + y's
        dft_bands_data = pd.DataFrame(
            np.hstack([self.dft_axis, self.dft_eigs]),
            columns=["x"] + self.dft_data_colnames
        )

        # append. 
        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": "dset_dftelbands",
                "data": dft_bands_data  # Store df as a single object in one row
            }
        ])
        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def add_wan_bands_data(self):
        # Read wannierqe eigs.
        wan_analysis = WannierQeAnalysis()
        wan_analysis.read_all()
        self.wan_eigs = wan_analysis.eigs - self.fermi_energy
        self.wan_axis = wan_analysis.wan_axis.reshape(-1, 1)
        self.wan_num_bands = self.wan_eigs.shape[1]

        # Create column names: y1, y2, ..., yN. 
        self.wan_data_colnames = [f"y{i+1}" for i in range(self.wan_num_bands)]


        # Build dataframe with x + y's
        wan_bands_data = pd.DataFrame(
            np.hstack([self.wan_axis, self.wan_eigs]),
            columns=["x"] + self.wan_data_colnames
        )

        # append.
        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": "dset_wannierqebands",
                "data": wan_bands_data  # Store df as a single object in one row
            }
        ])
        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def get_data(self):
        # Get kpath info.
        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.dft_axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Add the data.
        if os.path.exists('./dftelbands/dftelbands.xml'): self.add_dft_data()
        self.add_wan_bands_data()

    def add_dft_bands_figure(self):
        # Add the dft bands. 
        for ib in range(self.dft_num_bands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': f'wannier_bands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.dft_axis[0], self.dft_axis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Wannier Bandstructure',
                    'dset_name': 'dset_dftelbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.dft_data_colnames[ib]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': 'DFT' if ib==0 else None,
                },
            ])
            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

    def add_wannier_bands_figure(self):
        # Add the wannier bands.
        for ib in range(self.wan_num_bands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': f'wannier_bands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.dft_axis[0], self.dft_axis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Wannier Bandstructure',
                    'dset_name': 'dset_wannierqebands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.wan_data_colnames[ib]],
                    'color': 'red', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': 'Wannier' if ib==0 else None,
                },
            ])
            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

    def add_hr_figure(self):
        with h5py.File('./wan/wannier.h5', 'r') as hf:
            hr = hf['hr'][:]
            rpts = hf['hr_R'][:]

        nR = hr.shape[0]

        norm_list = []
        max_H_list = []
        for iR in range(nR):
            norm = np.linalg.norm(rpts[iR])
            norm_list.append(norm)

            max_H = np.max(np.abs(hr[iR].flatten()))
            max_H_list.append(max_H)

        norm_array = np.array(norm_list)
        max_H_array = np.array(max_H_list)

        # Sort norm_array and max_H_array based on norm_array
        sort_idx = np.argsort(norm_array)
        norm_array = norm_array[sort_idx]
        max_H_array = max_H_array[sort_idx]
        
        # Create df. 
        hr_df: pd.DataFrame = pd.DataFrame({
            'x': norm_array,
            'y': max_H_array,
        })

        # Add to self.dsets_df
        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": "dset_hr",
                "data": hr_df  # Store df as a single object in one row
            }
        ])
        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

        # Add scatter and line plot to figure. 
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': f'wannier_hr',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTER, 'axis': None,
                'xlabel': 'R (lattice units)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'max|H(R)| (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Wannier Hamiltonian Decay',
                'dset_name': 'dset_hr',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['y'],
                'color': 'blue', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
            {
                'fig_name': f'wannier_hr',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'R (lattice units)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'max|H(R)| (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Wannier Hamiltonian Decay',
                'dset_name': 'dset_hr',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['y'],
                'color': 'red', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
        ])
        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)

    def set_figures(self):
        if os.path.exists('./dftelbands/dftelbands.xml'): self.add_dft_bands_figure()
        self.add_wannier_bands_figure()
        self.add_hr_figure()

#endregion