
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

#endregion

#region variables
#endregion

#region functions
def pdos_label_from_filename(filename):
    pattern = r'.*?pdos.dat.pdos_atm#(?P<idx>\d+)\((?P<symbol>\w+)\)_wfc#\d+\((?P<orbital>.*)\)'

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
class PdosPlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)

        # Get pdos files. Filter if needed.
        self.pdos_files: list = []
        pattern = jmespath.search('pdos.plot.filter_regex', self.inputdict)
        for file in glob.glob('./pdos/pdos.dat.pdos_atm*'):
            if pattern is not None:
                if re.search(rf'{pattern}', file):
                    self.pdos_files.append(file)
            else:
                self.pdos_files.append(file)
        self.pdos_files = sorted(self.pdos_files)

        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get fermi energy
        tree = ET.parse('./scf/scf.xml')
        root = tree.getroot()
        fermi_energy = float(root.findall('.//fermi_energy')[0].text) * Hartree

        # Get name
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        for pdos_file in self.pdos_files:
            # Read file, skip header
            data = np.loadtxt(pdos_file, skiprows=1)
            energy = data[:, 0] - fermi_energy
            ldos = data[:, 1]
            dos_data = pd.DataFrame({'x': energy, 'y': ldos})
            dset_name = pdos_label_from_filename(pdos_file)
            append_dset_df = pd.DataFrame([
                {
                    "name": dset_name,
                    "data": dos_data,
                },
            ])
            self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        # Set colors.
        colors = plt.cm.tab20(np.linspace(0, 1, len(self.pdos_files)))

        # Individual figures
        for pdos_file, color in zip(self.pdos_files, colors):
            label = pdos_label_from_filename(pdos_file)
            
            append_fig_df = pd.DataFrame([
                {
                    'fig_name': f'pdos_{label}',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': 'Energy (eV)', 'xlim': (-10, 10), 'xticks': None, 'xtick_labels': None,
                    'ylabel': 'PDOS', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} PDOS ({label})',
                    'dset_name': label,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': ['y'],
                    'color': 'blue',
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])

            append_overlay_df = pd.DataFrame([
                {
                    'fig_name': f'pdos_overlay',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': 'Energy (eV)', 'xlim': (-10, 10), 'xticks': None, 'xtick_labels': None,
                    'ylabel': 'PDOS', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} PDOS',
                    'dset_name': label,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': ['y'],
                    'color': color,
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': label,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df, append_overlay_df], ignore_index=True)

#endregion