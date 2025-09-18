#region modules
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

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class PhbandsPlot(PlotBase):
    def __init__(
        self,
        infile='./struct.freq.gp',
        outfile_prefix='phbands_ph',
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.infilename: str = infile
        self.outfile_prefix: str = outfile_prefix
        
        self.get_data()
        self.set_figures()

    def get_data(self):
        data = np.loadtxt(self.infilename)
        self.phbands = data[:, 1:]
        self.numbands: int = self.phbands.shape[1]

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file().inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(self.phbands.shape[1])]

        # Build dataframe with x + y's
        df = pd.DataFrame(
            np.hstack([self.axis, self.phbands]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_phbands"],
            "data": [df]  # Store df as a single object in one row
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

    def set_figures(self):
        for ib in range(self.numbands):
            append_fig_df: pd.DataFrame = pd.DataFrame([{
                'fig_name': self.outfile_prefix,
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': None, 'xlim': None, 'xticks': self.xticks, 'xtick_labels': self.xtick_labels,
                'ylabel': r'Energy ($cm^{-1}$)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': f'{self.struct_name} Phonon Bandstructure',
                'dset_name': 'dset_phbands',
                'dset_axis_cols': 'x',        
                'dset_data_cols': [self.data_colnames[ib]],
                'color': 'blue', 
                'xgrid': True,
                'ygrid': False,
                'legend_label': None,
            }])

            self.figs_df = pd.concat([self.figs_df, append_fig_df], ignore_index=True)            

class PhonopyPlot(PhbandsPlot):
    def __init__(
        self,
        infile='./band.yaml',
        outfile_prefix='phbands_phonopy',
        **kwargs
    ):
        super().__init__(infile=infile, outfile_prefix=outfile_prefix, **kwargs)

    def get_data(self):
        with open(self.infilename) as f: data = yaml.safe_load(f)

        nk = len(data['phonon'])
        nb = len(data['phonon'][0]['band'])

        # fill phbands
        self.phbands = np.zeros(shape=(nk, nb), dtype='f8')
        for (k, b), value in np.ndenumerate(self.phbands):
            self.phbands[k, b] = data['phonon'][k]['band'][b]['frequency']*33.356        # Factor in cm^{-1}

        self.numbands: int = self.phbands.shape[1]
        print(self.phbands.shape)

        self.xaxis, self.xticks, self.xtick_labels = Kpath.from_yamlfile().even_spaced_axis
        self.axis = self.xaxis.reshape(-1, 1)

        # Get name.
        inputdict: dict = InputYaml.from_yaml_file().inputdict
        active_idx: int = jmespath.search('structures.active_idx', inputdict)
        self.struct_name: str = jmespath.search(f'structures.list[{active_idx}].name', inputdict)

        # Create column names: y1, y2, ..., yN
        self.data_colnames = [f"y{i+1}" for i in range(self.phbands.shape[1])]

        # Build dataframe with x + y's
        df = pd.DataFrame(
            np.hstack([self.axis, self.phbands]),
            columns=["x"] + self.data_colnames
        )

        append_dset_df: pd.DataFrame = pd.DataFrame({
            "name": ["dset_phbands"],
            "data": [df]  # Store df as a single object in one row
        })

        self.dsets_df = pd.concat([self.dsets_df, append_dset_df], ignore_index=True)

#endregion