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

from fpflow.structure.qe.qe_struct import QeStruct 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeConvergenceScfPlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.get_data_and_figures()
        self.add_kpt_bandgap_convergence()

    def add_kpt_bandgap_convergence(self):
        bandgap_df = pd.DataFrame({
            'ecut': pd.Series(dtype='Float64'),
            'kgrid': pd.Series(dtype='object'),
            'special_kpt': pd.Series(dtype='string'),
            'bandgap': pd.Series(dtype='Float64'),
        })

        dirnames = glob.glob('./convergence/qescf/dset_*')
        
        # If dirs to plot is specified, filter dirs.
        dirs_to_plot = jmespath.search('convergence.qescf.plot.dirs', self.inputdict)
        if dirs_to_plot is not None:
            dirnames = [d for d in dirnames if os.path.basename(d) in dirs_to_plot]

        dirnames.sort()

        current_dir = os.getcwd()

        for idir in dirnames:
            # Change to idir directory
            os.chdir(idir)

            # Load bandstructure data
            tree = ET.parse('./dftelbands.xml')
            root = tree.getroot()
            eig_nodes = root.findall('.//ks_energies/eigenvalues')
            fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
            num_kpts = len(eig_nodes)
            num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
            dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
            for kpt_idx, node in enumerate(eig_nodes):
                dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - fermi_energy

            dft_eigs = dft_eigs
            num_bands = dft_eigs.shape[1]

            
            # data = np.loadtxt('./bandstructure_inteqp.dat', skiprows=2)
            # num_bands = np.unique(data[:, 1]).size
            # eqp = data[:, 6].reshape(num_bands, -1).T

            # Get input.yaml info
            inputdict = InputYaml.from_yaml_file('./input.yaml').inputdict
            scf_ecut = float(jmespath.search('scf.ecut', inputdict))
            scf_kgrid = 'x'.join(map(str, jmespath.search('scf.kgrid', inputdict)))
            # Qestruct.
            max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
                xc=jmespath.search('scf.xc', self.inputdict),
                is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
            ))

            # Get axis info.
            self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
            self.axis = self.xaxis.reshape(-1, 1)

            # Get special k-point bandgap (first and last columns)
            for special_kpt_idx, kpt_label in zip(self.xticks, self.xtick_labels):
                kpt_index = np.argmin(np.abs(self.axis.flatten() - special_kpt_idx))
                vbm = dft_eigs[kpt_index, max_val_bands-1]
                cbm = dft_eigs[kpt_index, max_val_bands]
                bandgap = cbm - vbm

                # Append a single row in a way that preserves the existing dtypes
                bandgap_df.loc[len(bandgap_df)] = {
                    'ecut': scf_ecut,
                    'kgrid': scf_kgrid,
                    'special_kpt': kpt_label,
                    'bandgap': bandgap,
                }

            # Change back to original directory
            os.chdir(current_dir)

        # Add the datasets and figures.
        # ecut.
        # Use transform-based filtering instead of groupby.apply to avoid future deprecations
        bandgap_df['kgrid_val'] = bandgap_df['kgrid'].apply(lambda x: np.prod([int(i) for i in x.split('x')]))
        # Sort by ecut, kgrid_val (desc)
        bandgap_df_sorted = bandgap_df.sort_values(['ecut', 'kgrid_val'], ascending=[True, False])
        # For each ecut, select rows with max kgrid_val using transform
        max_kgrid_per_ecut = bandgap_df_sorted.groupby('ecut')['kgrid_val'].transform('max')
        ecut_df = (
            bandgap_df_sorted[bandgap_df_sorted['kgrid_val'] == max_kgrid_per_ecut]
            [['ecut', 'bandgap', 'special_kpt']]
            .reset_index(drop=True)
        )

        # kgrid.
        # For each kgrid group, select rows with max ecut
        bandgap_df_sorted = bandgap_df.sort_values(['kgrid_val', 'ecut'], ascending=[True, False])
        max_ecut_per_kgrid = bandgap_df_sorted.groupby('kgrid')['ecut'].transform('max')
        kgrid_df = (
            bandgap_df_sorted[bandgap_df_sorted['ecut'] == max_ecut_per_kgrid]
            [['kgrid', 'bandgap', 'special_kpt']]
            .reset_index(drop=True)
        )

        # Add ecut, and kgrid datasets to self.dsets_df and corresponding figures to self.figs_df.
        # Add ecut datasets and figures
        for special_kpt in ecut_df['special_kpt'].unique():
            df = ecut_df[ecut_df['special_kpt'] == special_kpt][['ecut', 'bandgap']].reset_index(drop=True)
            dset_name = f'ecut_{special_kpt}'
            append_dset_df = pd.DataFrame([{"name": dset_name, "data": df}])
            self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

            # Add corresponding figure
            fig_name = f'convqescf_ecut'
            append_fig_df = pd.DataFrame([{
            'fig_name': fig_name,
            'figure': None,
            'subplot_nrow': 1,
            'subplot_ncol': 1,
            'subplot_idx': 1,
            'plot_type': PlotType.LINE,
            'axis': None,
            'xlabel': 'ecut',
            'xlim': None,
            'xticks': None,
            'xtick_labels': None,
            'ylabel': 'Bandgap (eV)',
            'ylim': None,
            'yticks': None,
            'ytick_labels': None,
            'zlabel': None,
            'zlim': None,
            'zticks': None,
            'ztick_labels': None,
            'z_inc': None,
            'z_azim': None,
            'title': f'Bandgap vs ecut',
            'dset_name': dset_name,
            'dset_axis_cols': 'ecut',
            'dset_data_cols': ['bandgap'],
            'color': None,
            'linewidth': 1,
            'xgrid': True,
            'ygrid': True,
            'legend_label': None,
            }])
            append_scatter_df = pd.DataFrame([{
            'fig_name': fig_name,
            'figure': None,
            'subplot_nrow': 1,
            'subplot_ncol': 1,
            'subplot_idx': 1,
            'plot_type': PlotType.SCATTER,
            'axis': None,
            'xlabel': 'ecut',
            'xlim': None,
            'xticks': None,
            'xtick_labels': None,
            'ylabel': 'Bandgap (eV)',
            'ylim': None,
            'yticks': None,
            'ytick_labels': None,
            'zlabel': None,
            'zlim': None,
            'zticks': None,
            'ztick_labels': None,
            'z_inc': None,
            'z_azim': None,
            'title': f'Bandgap vs ecut',
            'dset_name': dset_name,
            'dset_axis_cols': 'ecut',
            'dset_data_cols': ['bandgap'],
            'color': None,
            'linewidth': 1,
            'xgrid': True,
            'ygrid': True,
            'legend_label': f'{special_kpt}',
            }])
            self.figs_df = pd.concat([self.figs_df, append_fig_df, append_scatter_df], ignore_index=True)

        # Add kgrid datasets and figures
        for special_kpt in kgrid_df['special_kpt'].unique():
            df = kgrid_df[kgrid_df['special_kpt'] == special_kpt][['kgrid', 'bandgap']].reset_index(drop=True)
            dset_name = f'kgrid_{special_kpt}'
            append_dset_df = pd.DataFrame([{"name": dset_name, "data": df}])
            self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

            # Add corresponding figure
            fig_name = f'convqescf_kgrid'
            append_fig_df = pd.DataFrame([{
            'fig_name': fig_name,
            'figure': None,
            'subplot_nrow': 1,
            'subplot_ncol': 1,
            'subplot_idx': 1,
            'plot_type': PlotType.LINE,
            'axis': None,
            'xlabel': 'kgrid',
            'xlim': None,
            'xticks': None,
            'xtick_labels': None,
            'ylabel': 'Bandgap (eV)',
            'ylim': None,
            'yticks': None,
            'ytick_labels': None,
            'zlabel': None,
            'zlim': None,
            'zticks': None,
            'ztick_labels': None,
            'z_inc': None,
            'z_azim': None,
            'title': f'Bandgap vs kgrid',
            'dset_name': dset_name,
            'dset_axis_cols': 'kgrid',
            'dset_data_cols': ['bandgap'],
            'color': None,
            'linewidth': 1,
            'xgrid': True,
            'ygrid': True,
            'legend_label': None,
            }])
            append_scatter_df = pd.DataFrame([{
            'fig_name': fig_name,
            'figure': None,
            'subplot_nrow': 1,
            'subplot_ncol': 1,
            'subplot_idx': 1,
            'plot_type': PlotType.SCATTER,
            'axis': None,
            'xlabel': 'kgrid',
            'xlim': None,
            'xticks': None,
            'xtick_labels': None,
            'ylabel': 'Bandgap (eV)',
            'ylim': None,
            'yticks': None,
            'ytick_labels': None,
            'zlabel': None,
            'zlim': None,
            'zticks': None,
            'ztick_labels': None,
            'z_inc': None,
            'z_azim': None,
            'title': f'Bandgap vs kgrid',
            'dset_name': dset_name,
            'dset_axis_cols': 'kgrid',
            'dset_data_cols': ['bandgap'],
            'color': None,
            'linewidth': 1,
            'xgrid': True,
            'ygrid': True,
            'legend_label': f'{special_kpt}',
            }])
            self.figs_df = pd.concat([self.figs_df, append_fig_df, append_scatter_df], ignore_index=True)
     
    def get_from_single_folder(self, dest_dir: str, color, src_dir: str=os.getcwd()):
        # Change to dest dir.
        os.chdir(dest_dir)

        tree = ET.parse('./dftelbands.xml')
        root = tree.getroot()

        eig_nodes = root.findall('.//ks_energies/eigenvalues')
        fermi_energy = float(root.findall('.//fermi_energy')[0].text)*Hartree
        num_kpts = len(eig_nodes)
        num_bands = np.fromstring(eig_nodes[0].text, sep=' ', dtype='f8').size
        dft_eigs = np.zeros(shape=(num_kpts, num_bands), dtype='f8')
        for kpt_idx, node in enumerate(eig_nodes):
            dft_eigs[kpt_idx, :] = np.fromstring(node.text, sep=' ', dtype='f8')*Hartree - fermi_energy

        self.dft_eigs = dft_eigs
        self.numbands = dft_eigs.shape[1]

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file('./input.yaml').inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)
        scf_kgrid: str = 'x'.join(list(map(str, jmespath.search('scf.kgrid', inputdict))))
        scf_ecut: str = 'ecut' + str(jmespath.search('scf.ecut', inputdict))
        dset_name: str = f'dset_{scf_ecut}_{scf_kgrid}'

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(self.dft_eigs.shape[1])]

        # Build dataframe with x + y's
        df = pd.DataFrame(
            np.hstack([self.axis, self.dft_eigs]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": [dset_name],
            "data": [df]  # Store df as a single object in one row
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

        for ib in range(self.numbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([{
                'fig_name': 'convqescf',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': 'Energy (eV)', 'ylim': (-10, 10), 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} DFT Bandstructure',
                'dset_name': dset_name,
                'dset_axis_cols': 'x',        
                'dset_data_cols': [self.data_colnames[ib]],
                'color': color, 
                'linewidth': 1, 
                'xgrid': True,
                'ygrid': False,
                'legend_label': f'{scf_ecut}_{scf_kgrid}' if ib==0 else None,
            }])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)    

        # Change back to src dir.
        os.chdir(src_dir)
                  
    def get_data_and_figures(self):
        dirs = glob.glob('./convergence/qescf/dset_*')

        # If dirs to plot is specified, filter dirs.
        dirs_to_plot = jmespath.search('convergence.qescf.plot.dirs', self.inputdict)
        if dirs_to_plot is not None:
            dirs = [d for d in dirs if os.path.basename(d) in dirs_to_plot]

        dirs.sort()

        colors = plt.cm.tab20(np.linspace(0, 1, len(dirs)))

        for idir, color in zip(dirs, colors):
            self.get_from_single_folder(idir, color)

#endregion