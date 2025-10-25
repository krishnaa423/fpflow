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

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwZdPlot(PlotBase):
    def __init__(
        self,
        **kwargs,
    ):
        super().__init__(**kwargs)
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        # Get name
        inputdict: dict = self.inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Get kpts.
        Qpts: list = Kpts.from_kgrid(
            kgrid=[
                jmespath.search('bseq.qgrid[0]', self.inputdict),
                jmespath.search('bseq.qgrid[1]', self.inputdict),
                jmespath.search('bseq.qgrid[2]', self.inputdict),
            ],
            is_reduced=False,
        ).bseq_qpts

        # Get kpath.
        self.kpath: Kpath = Kpath.from_yamlfile()
        self.xaxis, self.xticks, self.xtick_labels = self.kpath.even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get phgrid.

        # Get phbands. 
        data = np.loadtxt('./struct.freq.gp')
        self.phbands = data[:, 1:]
        self.phbands *= 0.123984   # Convert to meV. 1 cm-1 = 0.123984 meV
        self.numbands: int = self.phbands.shape[1]


        # Get exgrid and exbands.
        # exband.fmt: rows = Qpts, cols = bands
        import scipy.interpolate
        exband = np.loadtxt('./zd/epw/G_full_epmatq/exband.fmt')
        exgrid = np.array(Qpts)  # shape (n_qpts, 3)
        kpath_kpts = np.array(self.kpath.kpts)  # shape (n_kpts, 3)
        num_bands = exband.shape[1]
        exbands = np.zeros((kpath_kpts.shape[0], num_bands))
        # Use RBFInterpolator for robust interpolation and extrapolation
        for ib in range(num_bands):
            rbf_interp = scipy.interpolate.RBFInterpolator(exgrid, exband[:, ib], neighbors=10)
            exbands[:, ib] = rbf_interp(kpath_kpts)
        self.exgrid = exgrid
        self.exbands = exbands
        self.num_exbands = num_bands

        # Get xctph.

        # Get asq.
        # Read Asq.plrn and interpolate onto exciton bands
        asq_path = './zd/epw/Asq.plrn'
        asq_coords = []  # (Qpt, band_idx)
        asq_real = []
        asq_imag = []
        with open(asq_path, 'r') as f:
            lines = f.readlines()
        # Skip header (first line)
        for line in lines[1:]:
            parts = line.split()
            if len(parts) == 4:
                grid_idx = int(parts[0]) - 1  # 1-based to 0-based
                band_idx = int(parts[1]) - 1
                real_val = float(parts[2])
                imag_val = float(parts[3])
                # Use Qpts[grid_idx] for actual k-point coordinates
                asq_coords.append(list(Qpts[grid_idx]) + [band_idx])
                asq_real.append(real_val)
                asq_imag.append(imag_val)
        asq_coords = np.array(asq_coords)
        asq_real = np.array(asq_real)
        asq_imag = np.array(asq_imag)

        # Interpolate onto (kpath_kpt, band_idx)
        n_kpts = self.exbands.shape[0]
        n_bands = self.exbands.shape[1]
        interp_points = np.array([list(kpath_kpts[i]) + [b] for i in range(n_kpts) for b in range(n_bands)])
        rbf_real = scipy.interpolate.RBFInterpolator(asq_coords, asq_real, neighbors=10)
        rbf_imag = scipy.interpolate.RBFInterpolator(asq_coords, asq_imag, neighbors=10)
        asq_real_interp = rbf_real(interp_points).reshape(n_kpts, n_bands)
        asq_imag_interp = rbf_imag(interp_points).reshape(n_kpts, n_bands)
        self.asq_real = asq_real_interp
        self.asq_imag = asq_imag_interp
        self.asq_mag = np.sqrt(asq_real_interp**2 + asq_imag_interp**2)

        # Get bmat.
        # Read bmat.band.plrn and interpolate weights onto kpath kpts
        bmat_path = './zd/epw/bmat.band.plrn'
        bmat_coords = []  # (Qpt, phonon_band)
        bmat_weights = []
        with open(bmat_path, 'r') as f:
            bmat_lines = f.readlines()
        # Skip header (first line)
        for line in bmat_lines[1:]:
            parts = line.split()
            if len(parts) >= 3:
                grid_idx = int(parts[0]) - 1  # 1-based to 0-based
                ph_band = int(parts[1]) - 1
                weight = float(parts[-1])
                bmat_coords.append(list(Qpts[grid_idx]) + [ph_band])
                bmat_weights.append(weight)
        bmat_coords = np.array(bmat_coords)
        bmat_weights = np.array(bmat_weights)

        # Interpolate onto (kpath_kpt, phonon_band)
        n_kpts_bmat = kpath_kpts.shape[0]
        n_phbands = np.max([int(line.split()[1]) for line in bmat_lines[1:] if len(line.split()) >= 2])
        interp_points_bmat = np.array([list(kpath_kpts[i]) + [b] for i in range(n_kpts_bmat) for b in range(n_phbands)])
        rbf_bmat = scipy.interpolate.RBFInterpolator(bmat_coords, bmat_weights, neighbors=10)
        bmat_weights_interp = rbf_bmat(interp_points_bmat).reshape(n_kpts_bmat, n_phbands)
        self.bmat_weights = bmat_weights_interp

        # Get asq on exbands.

        # Get bmat on phbands.

        # exbands df 
        self.exband_colnames = [f"y{i+1}" for i in range(self.exbands.shape[1])]
        self.exband_weight_colnames = [f"s{i+1}" for i in range(self.exbands.shape[1])]
        df = pd.DataFrame(
            np.hstack([self.axis, self.exbands, self.asq_mag*100]),
            columns=["x"] + self.exband_colnames + self.exband_weight_colnames
        )
        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_exbands"],
            "data": [df]  # Store df as a single object in one row
        })
        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)
        df.to_csv('debug_epw_exbands.csv', index=False)

        # --- phbands df ---
        self.phband_colnames = [f"y{i+1}" for i in range(self.phbands.shape[1])]
        self.phband_weight_colnames = [f"s{i+1}" for i in range(self.phbands.shape[1])]
        df_ph = pd.DataFrame(
            np.hstack([self.axis, self.phbands, self.bmat_weights*30]),
            columns=["x"] + self.phband_colnames + self.phband_weight_colnames
        )
        append_ph_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_phbands"],
            "data": [df_ph]
        })
        self.dsets_df = pd.concat([self.dsets_df, append_ph_dset_df], ignore_index=True)

    def set_figures(self):
        for ib in range(self.num_exbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': 'zd_exbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} Exciton Bandstructure',
                    'dset_name': 'dset_exbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.exband_colnames[ib]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'zd_asq',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Projection on Exciton Bands',
                    'dset_name': 'dset_exbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.exband_colnames[ib]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'zd_asq',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Energy (eV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Projection on Exciton Bands',
                    'dset_name': 'dset_exbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.exband_colnames[ib], self.exband_weight_colnames[ib]],
                    'color': 'red', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)            

        # Add figures for phbands
        for ib in range(self.phbands.shape[1]):
            append_ph_fig_df: pd.DataFrame = pd.DataFrame([
                {
                    'fig_name': 'zd_phbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.LINE, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Phonon Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Displacement projected on Phonon Bandstructure',
                    'dset_name': 'dset_phbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.phband_colnames[ib]],
                    'color': 'blue', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
                {
                    'fig_name': 'zd_phbands',
                    'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                    'plot_type': PlotType.SCATTER, 'axis': None,
                    'xlabel': None, 'xlim': (self.xaxis[0], self.xaxis[-1]), 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                    'ylabel': 'Phonon Energy (meV)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                    'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                    'z_inc': None, 'z_azim': None,
                    'title': f'{self.struct_name} STE Displacement projected on Phonon Bandstructure',
                    'dset_name': 'dset_phbands',
                    'dset_axis_cols': 'x',        
                    'dset_data_cols': [self.phband_colnames[ib], self.phband_weight_colnames[ib]],
                    'color': 'red', 
                    'xgrid': True,
                    'ygrid': False,
                    'legend_label': None,
                },
            ])
            self.figs_df = pd.concat([self.figs_df, append_ph_fig_df], ignore_index=True)

#endregion