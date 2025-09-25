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
        self.add_eig_convergence_figure()

    def _add_by_col(self, colname: str):
        if not self.eig_df.empty:
            # Pivot the dataframe so each eig_idx becomes its own column
            pivot_df = self.eig_df.pivot_table(
            index=colname,
            columns='eig_idx',
            values='eig_val',
            aggfunc='max'
            ).reset_index()

            # Rename columns: eig_{eig_idx}
            pivot_df = pivot_df.rename(
            columns={col: f"eig_{col}" for col in pivot_df.columns if isinstance(col, (int, np.integer))}
            )
            # x column is coarse_kgrid
            pivot_df = pivot_df.rename(columns={colname: 'x'})

            # Add the dataset to dsets_df
            dset_name = f'eig_convergence_{colname}'
            append_eig_conv_df = pd.DataFrame({
            "name": [dset_name],
            "data": [pivot_df]
            })
            self.dsets_df = pd.concat([self.dsets_df, append_eig_conv_df], ignore_index=True)

            eig_unique = self.eig_df['eig_idx'].unique()

            # Add a figure for each eig_idx
            for eig_idx in eig_unique:
                eig_col = f"eig_{eig_idx}"
                append_fig_df = pd.DataFrame([{
                    'fig_name': f'convbgwbse_eig_{colname}',
                    'figure': None,
                    'subplot_nrow': 1,
                    'subplot_ncol': 1,
                    'subplot_idx': 1,
                    'plot_type': PlotType.LINE,
                    'axis': None,
                    'xlabel': f'{colname}',
                    'xlim': None,
                    'xticks': None,
                    'xtick_labels': None,
                    'ylabel': 'Eigenvalue (eV)',
                    'ylim': None,
                    'yticks': None,
                    'ytick_labels': None,
                    'zlabel': None,
                    'zlim': None,
                    'zticks': None,
                    'ztick_labels': None,
                    'z_inc': None,
                    'z_azim': None,
                    'title': f'BSE Eigenvalue Convergence',
                    'dset_name': dset_name,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': [eig_col],
                    'color': None,
                    'linewidth': 2,
                    'xgrid': True,
                    'ygrid': True,
                    'legend_label': f'eig_{eig_idx}',
                }])

                append_scatter_df = pd.DataFrame([{
                    'fig_name': f'convbgwbse_eig_{colname}',
                    'figure': None,
                    'subplot_nrow': 1,
                    'subplot_ncol': 1,
                    'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER,
                    'axis': None,
                    'xlabel': f'{colname}',
                    'xlim': None,
                    'xticks': None,
                    'xtick_labels': None,
                    'ylabel': 'Eigenvalue (eV)',
                    'ylim': None,
                    'yticks': None,
                    'ytick_labels': None,
                    'zlabel': None,
                    'zlim': None,
                    'zticks': None,
                    'ztick_labels': None,
                    'z_inc': None,
                    'z_azim': None,
                    'title': f'BSE Eigenvalue Convergence',
                    'dset_name': dset_name,
                    'dset_axis_cols': 'x',
                    'dset_data_cols': [eig_col],
                    'color': None,
                    'linewidth': 2,
                    'xgrid': True,
                    'ygrid': True,
                    'legend_label': None,
                }])

                self.figs_df = pd.concat([self.figs_df, append_scatter_df, append_fig_df], ignore_index=True)

    def add_eig_convergence_figure(self):
        '''
        1). Loop over folders. create the df with coarse grid, fine grid, c, v, eig_idx, eig_val entry for each. 
        2). x axis: coarse grid, yaxis: eig val, hue: eig idx. 
        '''
        
        dirs = glob.glob('./convergence/bgwbse/dset_*')

        # If dirs to plot is specified, filter dirs.
        dirs_to_plot = jmespath.search('convergence.bgwbse.plot.dirs', self.inputdict)
        if dirs_to_plot is not None:
            dirs = [d for d in dirs if os.path.basename(d) in dirs_to_plot]

        dirs.sort()
        current_dir = os.getcwd()

        # 1). Get the data.
        self.eig_df = pd.DataFrame({
            'bse_nv': pd.Series(dtype='int'),
            'bse_nc': pd.Series(dtype='int'),
            'coarse_kgrid': pd.Series(dtype='str'),
            'fine_kgrid': pd.Series(dtype='str'),
            'eig_idx': pd.Series(dtype='Int32'),
            'eig_val': pd.Series(dtype='Float64'),
        })

        for idir in dirs:
            os.chdir(idir)
            inputdict: dict = InputYaml.from_yaml_file('./input.yaml').inputdict
            active_idx: int = jmespath.search('structures.active_idx', inputdict)
            struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)
            bse_nv: int = jmespath.search('bse.absorption.val_bands', inputdict)
            bse_nc: int = jmespath.search('bse.absorption.cond_bands', inputdict)
            bse_coarse_kgrid: str = 'x'.join(list(map(str, jmespath.search('wfn.kgrid', inputdict))))
            bse_fine_grid: str = 'x'.join(list(map(str, jmespath.search('wfnfi.kgrid', inputdict))))

            eigfile_data = np.loadtxt('./eigenvalues.dat', dtype='f8', skiprows=4)
            eigs = eigfile_data[:, 0]  # in eV

            # Get eig idx to plot.
            eigs_to_plot: list = jmespath.search('convergence.bgwbse.plot.eigs', self.inputdict)
            if eigs_to_plot is None:
                eigs_to_plot = [0]

            for eig_idx in eigs_to_plot:
                if eig_idx < 0 or eig_idx >= len(eigs):
                    raise ValueError(f'Error: Specified eig index {eig_idx} is out of range [0, {len(eigs)-1}]')

                append_eig_df: pd.DataFrame = pd.DataFrame({
                    "bse_nv": [bse_nv],
                    "bse_nc": [bse_nc],
                    "coarse_kgrid": [bse_coarse_kgrid],
                    "fine_kgrid": [bse_fine_grid],
                    "eig_idx": [eig_idx],
                    "eig_val": [eigs[eig_idx]],
                })
                
                frames = [self.eig_df, append_eig_df]
                # Keep only DataFrames that have at least one non-NA value
                frames = [df for df in frames if not df.empty and not df.isna().all().all()]

                if frames:
                    self.eig_df = pd.concat(frames, ignore_index=True)

            os.chdir(current_dir)

        # 2). Create the figure. 
        self._add_by_col('coarse_kgrid')
        self._add_by_col('fine_kgrid')
        
    def get_from_single_folder(self, dest_dir: str, color, src_dir: str=os.getcwd()):
        # Change to dest dir.
        os.chdir(dest_dir)

        abs_eh_data = np.loadtxt('./absorption_eh.dat', dtype='f8', skiprows=4)
        axis = abs_eh_data[:, 0]
        eh_data = abs_eh_data[:, 1]

        # Load eigenvalue data.
        eigfile_data = np.loadtxt('./eigenvalues.dat', dtype='f8', skiprows=4)
        eigs = eigfile_data[:, 0]
        dipole = eigfile_data[:, 1]

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file('./input.yaml').inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)
        bse_nv: str = 'v' + str(jmespath.search('bse.absorption.val_bands', inputdict))
        bse_nc: str = 'c' + str(jmespath.search('bse.absorption.cond_bands', inputdict))
        bse_coarse_kgrid: str = 'c' + 'x'.join(list(map(str, jmespath.search('wfn.kgrid', inputdict))))
        bse_fine_grid: str = 'f' + 'x'.join(list(map(str, jmespath.search('wfnfi.kgrid', inputdict))))
        dset_name: str = f'dset_{bse_nv}_{bse_nc}_{bse_coarse_kgrid}_{bse_fine_grid}'

        append_abs_df: pd.DataFrame = pd.DataFrame({
            "name": [dset_name],
            "data": [pd.DataFrame({
                "x": axis,
                "y_eh": eh_data,
            })]
        })

        append_dipole_df: pd.DataFrame = pd.DataFrame({
            "name": [dset_name + '_dipole'],
            "data": [pd.DataFrame({
                "x": eigs,
                "y_dipole": dipole,
            })]
        })

        self.dsets_df = pd.concat([self.dsets_df, append_abs_df, append_dipole_df], ignore_index=True)

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
                'legend_label': f'{bse_nv}_{bse_nc}_{bse_coarse_kgrid}_{bse_fine_grid}',
            },
        ])

        append_dipole_fig_df: pd.DataFrame = pd.DataFrame([
            {
                'fig_name': 'convqebse_dipole',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.STEM, 'axis': None,
                'xlabel': 'Energy (eV)', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': r'Dipole Matrix Element (arb.)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} BSE Dipole Matrix Elements',
                'dset_name': dset_name + '_dipole',
                'dset_axis_cols': 'x',        
                'dset_data_cols': ['y_dipole'],
                'color': color, 
                'linewidth': 1, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': f'{bse_nv}_{bse_nc}_{bse_coarse_kgrid}_{bse_fine_grid}',
            },
        ])

        self.figs_df = pd.concat([self.figs_df, append_fig_df, append_dipole_fig_df], ignore_index=True)   

        # Change back to src dir.
        os.chdir(src_dir)

    def get_data_and_figures(self):
        dirs = glob.glob('./convergence/bgwbse/dset_*')
        
        # If dirs to plot is specified, filter dirs.
        dirs_to_plot = jmespath.search('convergence.bgwbse.plot.dirs', self.inputdict)
        if dirs_to_plot is not None:
            dirs = [d for d in dirs if os.path.basename(d) in dirs_to_plot]

        dirs.sort()

        colors = plt.cm.tab20(np.linspace(0, 1, len(dirs)))

        for idir, color in zip(dirs, colors):
            self.get_from_single_folder(idir, color)

#endregion