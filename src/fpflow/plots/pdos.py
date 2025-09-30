#region modules
from fpflow.plots.plot import PlotBase, PlotType
import pandas as pd
from ase.units import Hartree, eV
import jmespath
import xml.etree.ElementTree as ET
import numpy as np 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class PdosPlot(PlotBase):
    def __init__(
        self,
        infile_scf='./scf.xml',
        infile_dos='./struct_pdos.dat',
        outfile_prefix='pdos',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.infile_scf: str = infile_scf
        self.infile_dos: str = infile_dos
        self.outfile_prefix: str = outfile_prefix
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get fermi energy. 
        tree = ET.parse(self.infile_scf)
        root = tree.getroot()
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree

        data = np.loadtxt(self.infile_dos, skiprows=1)
        self.axis = (data[:, 0] - fermi_energy).reshape(-1, 1)
        self.dos = data[:, 1].reshape(-1, 1)

        # Get name.
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Build dataframe with x + y's
        dos_data = pd.DataFrame(
            np.concat([self.axis, self.dos], axis=1),
            columns=["x"] + ['y']
        )

        append_dset_df: pd.DataFrame = pd.DataFrame([
            {
                "name": "dset_dos",
                "data": dos_data  # Store df as a single object in one row
            },
        ])

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': self.outfile_prefix,
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'Energy (eV)', 'xlim': (-10, 10), 'xticks': None, 'xtick_labels': None,
                'ylabel': 'DOS', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} DOS',
                'dset_name': 'dset_dos',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['y'],
                'color': 'blue', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            },
        ])

        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)  

#endregion