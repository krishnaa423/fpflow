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
import glob 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwConvergenceBsePlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.get_data_and_figures()

    def get_from_single_folder(self, dest_dir: str, color, src_dir: str=os.getcwd()):
        # Change to dest dir.
        os.chdir(dest_dir)

        abs_eh_data = np.loadtxt('./absorption_eh.dat', dtype='f8', skiprows=4)
        axis = abs_eh_data[:, 0]
        eh_data = abs_eh_data[:, 1]

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file('./input.yaml').inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)
        bse_nv: str = 'v' + str(jmespath.search('bse.absorption.val_bands', inputdict))
        bse_nc: str = 'c' + str(jmespath.search('bse.absorption.cond_bands', inputdict))
        bse_kgrid: str = 'x'.join(list(map(str, jmespath.search('wfn.kgrid', inputdict))))
        dset_name: str = f'dset_{bse_nv}_{bse_nc}_{bse_kgrid}'

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": [dset_name],
            "data": [pd.DataFrame({
                "x": axis,
                "y_eh": eh_data,
            })]
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

        append_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'convqebse',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'Energy (eV)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': r'$\epsilon_2(\omega) (arb.)$', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} BSE Absorption',
                'dset_name': dset_name,
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['y_eh'],
                'color': color, 
                'linewidth': 1, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': f'{bse_nv}_{bse_nc}_{bse_kgrid}',
            },
        ])

        self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)   


        # Change back to src dir.
        os.chdir(src_dir)

    def get_data_and_figures(self):
        dirs = glob.glob('./convergence/qe/bse/*')
        dirs.sort()

        colors = plt.cm.tab20(np.linspace(0, 1, len(dirs)))

        for idir, color in zip(dirs, colors):
            self.get_from_single_folder(idir, color)

#endregion