
#region modules
from fpflow.plots.plot import PlotBase, PlotType
import pandas as pd
from ase.units import Hartree, eV
import jmespath
import xml.etree.ElementTree as ET
import numpy as np
import glob
import os
import re
import matplotlib.pyplot as plt 
from fpflow.structure.kpath import Kpath

#endregion

#region variables
#endregion

#region functions
def kpdos_label_from_filename(filename):
    pattern = r'.*?struct_kpdos.dat.pdos_atm#(?P<idx>\d+)\((?P<symbol>\w+)\)_wfc#\d+\((?P<orbital>.*)\)'

    match = re.match(pattern, filename)
    if match:
        atom_idx = match.group('idx')
        atom_symbol = match.group('symbol')
        orbital = match.group('orbital')
        label = f'{atom_idx}{atom_symbol}_{orbital}'
        return label
    else:
        raise ValueError(f"Filename '{filename}' does not match the expected pattern.")

#endregion

#region classes
class KpdosPlot(PlotBase):
    def __init__(
        self,
        infile_scf='./scf.xml',
        kpdos_glob='./struct_kpdos.dat.pdos_atm*',
        outfile_prefix='kpdos',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.infile_scf: str = infile_scf
        self.kpdos_glob: str = kpdos_glob
        self.outfile_prefix: str = outfile_prefix
        self.kpdos_files = sorted(glob.glob(self.kpdos_glob))
        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get fermi energy
        tree = ET.parse(self.infile_scf)
        root = tree.getroot()
        fermi_energy = float(root.findall('.//fermi_energy')[0].text) * Hartree

        # Get name
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        scale_factor: float = jmespath.search('kpdos.plot.scatter_scale_factor', inputdict)
        if scale_factor is None:
            scale_factor = 1.0

        for kpdos_file in self.kpdos_files:
            # Read file, skip header
            data = np.loadtxt(kpdos_file, skiprows=1)
            kpts = (data[:, 0] - 1).astype(float)  # Convert to zero-based index
            energy = data[:, 1] - fermi_energy
            ldos = data[:, 2] * scale_factor
            kpdos_data = pd.DataFrame({'x': kpts, 'y': energy, 'size': ldos})
            dset_name = kpdos_label_from_filename(kpdos_file)
            append_dset_df = pd.DataFrame([
                {
                    "name": dset_name,
                    "data": kpdos_data,
                },
            ])
            self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        # Set colors.
        colors = plt.cm.tab20(np.linspace(0, 1, len(self.kpdos_files)))

        # Individual figures
        for kpdos_file, color in zip(self.kpdos_files, colors):
            label = kpdos_label_from_filename(kpdos_file)

            append_fig_df = pd.DataFrame([
                {
                    'fig_name': f'kpdos_{label}',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} k-resolved PDOS ({label})',
                    'dset_name': label,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': ['y', 'size'],
                    'color': 'blue',
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])

            append_overlay_df = pd.DataFrame([
                {
                    'fig_name': f'kpdos_overlay',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} k-resolved PDOS',
                    'dset_name': label,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': ['y', 'size'],
                    'color': color,
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': label,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df, append_overlay_df], ignore_index=True)

#endregion