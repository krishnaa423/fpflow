#region modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
from matplotlib.axis import Axis 
from matplotlib.figure import Figure
from cycler import cycler
import pandas as pd 
import h5py
import os 

#endregion

#region variables
aps_rcparams = {
    # --------- Fonts & text ---------
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.5,              # base size ~9pt
    "axes.titlesize": 9.5,
    "axes.labelsize": 9.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 8.0,
    "figure.titlesize": 10.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    # --------- Figure layout ---------
    "figure.figsize": (3.4, 2.3),  # APS single col ~3.4in wide
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.format": "pdf",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.01,
    "savefig.transparent": False,
    "figure.constrained_layout.use": True,

    # --------- Lines, markers, errorbars ---------
    "lines.linewidth": 1.0,
    "lines.markersize": 3.0,
    "lines.markeredgewidth": 0.5,
    "lines.antialiased": True,
    "errorbar.capsize": 2.0,

    # --------- Patches / fills ---------
    "patch.linewidth": 0.8,
    "hatch.linewidth": 0.4,

    # --------- Axes & spines ---------
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.labelpad": 2.5,
    "axes.titlepad": 4.0,
    "axes.grid": False,
    "axes.axisbelow": True,
    "axes.formatter.use_mathtext": True,
    "axes.formatter.limits": (-3, 3),
    "axes.spines.top": True,
    "axes.spines.right": True,

    # --------- Ticks ---------
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.0,
    "xtick.minor.size": 1.5,
    "xtick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "xtick.minor.visible": True,
    "ytick.major.size": 3.0,
    "ytick.minor.size": 1.5,
    "ytick.major.width": 0.8,
    "ytick.minor.width": 0.6,
    "ytick.minor.visible": True,

    # --------- Legends ---------
    "legend.frameon": False,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.4,
    "legend.borderaxespad": 0.5,
    "legend.labelspacing": 0.3,
    "legend.borderpad": 0.3,

    # --------- Images / fields ---------
    "image.cmap": "viridis",
    "image.interpolation": "nearest",
    "image.origin": "lower",
    "contour.negative_linestyle": "dashed",
    # "contour.linewidths": 0.8,

    # --------- Boxplots ---------
    "boxplot.vertical": True,
    "boxplot.whiskerprops.linewidth": 0.8,
    "boxplot.boxprops.linewidth": 0.8,
    "boxplot.capprops.linewidth": 0.8,
    "boxplot.medianprops.linewidth": 0.9,
    "boxplot.flierprops.marker": "o",
    "boxplot.flierprops.markersize": 2.5,

    # --------- Grids (off by default) ---------
    "grid.linewidth": 0.6,
    "grid.alpha": 0.3,
    "grid.linestyle": "-",
    "grid.color": "#b0b0b0",

    # --------- Performance ---------
    "path.simplify": True,
    "path.simplify_threshold": 0.0,
    "agg.path.chunksize": 20000,

    # --------- Interactivity ---------
    "figure.autolayout": False,
}
nature_rcparams = {
    # --------- Fonts & text ---------
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 8.0,
    "figure.titlesize": 10.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    # --------- Figure layout ---------
    "figure.figsize": (85/25.4, (85/25.4) * 0.68),  # ≈ 3.35 x 2.28 in
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.format": "pdf",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.transparent": False,
    "figure.constrained_layout.use": True,

    # --------- Lines, markers, errorbars ---------
    "lines.linewidth": 1.1,
    "lines.markersize": 3.8,
    "lines.markeredgewidth": 0.6,
    "lines.dash_capstyle": "butt",
    "lines.solid_capstyle": "butt",
    "lines.antialiased": True,
    "errorbar.capsize": 2.5,
    "patch.linewidth": 0.8,
    "patch.antialiased": True,
    "hatch.linewidth": 0.4,

    # --------- Axes & spines ---------
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.labelpad": 2.5,
    "axes.titlepad": 4.0,
    "axes.grid": False,
    "axes.axisbelow": True,
    "axes.prop_cycle": (
        cycler("color", [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
        ])
    ),
    "axes.formatter.use_mathtext": True,
    "axes.formatter.limits": (-3, 3),
    "axes.spines.top": True,
    "axes.spines.right": True,

    # --------- Ticks ---------
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.2,
    "xtick.minor.size": 1.8,
    "xtick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "xtick.minor.visible": True,
    "ytick.major.size": 3.2,
    "ytick.minor.size": 1.8,
    "ytick.major.width": 0.8,
    "ytick.minor.width": 0.6,
    "ytick.minor.visible": True,

    # --------- Legends ---------
    "legend.frameon": False,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.4,
    "legend.borderaxespad": 0.5,
    "legend.labelspacing": 0.3,
    "legend.borderpad": 0.3,

    # --------- Images / fields ---------
    "image.cmap": "viridis",
    "image.interpolation": "nearest",
    "image.origin": "lower",
    "contour.negative_linestyle": "dashed",
    "contour.linewidth": 0.8,

    # --------- Distributions ---------
    "hist.bins": 50,
    "boxplot.vertical": True,
    "boxplot.whiskerprops.linewidth": 0.8,
    "boxplot.boxprops.linewidth": 0.8,
    "boxplot.capprops.linewidth": 0.8,
    "boxplot.medianprops.linewidth": 0.9,
    "boxplot.flierprops.marker": "o",
    "boxplot.flierprops.markersize": 2.8,
    "boxplot.flierprops.markeredgewidth": 0.4,

    # --------- 3D ---------
    "grid.linewidth": 0.6,
    "grid.alpha": 0.3,

    # --------- Grids ---------
    "grid.linestyle": "-",
    "grid.color": "#b0b0b0",

    # --------- Performance & vector safety ---------
    "path.simplify": True,
    "path.simplify_threshold": 0.0,
    "agg.path.chunksize": 20000,

    # --------- Interactivity ---------
    "figure.autolayout": False,
}

#endregion

#region functions
def set_common_rcparams(theme='nature'):
    match theme:
        case 'aps':
            mpl.rcParams.update(aps_rcparams)
        case 'nature':
            mpl.rcParams.update(nature_rcparams)
        case _:
            plt.style.use('seaborn-v0_8-whitegrid')

#endregion

#region classes
class AxisPlot:
    def __init__(
        self,
        axis: Axis,
        dset_row: pd.Series,
        figs_row: pd.Series,
        **kwargs,
    ):
        self.axis: Axis = axis
        self.dset_row: pd.Series = dset_row
        self.figs_row: pd.Series = figs_row

        for key, value in kwargs.items():
            setattr(self, key, value)

    def plot(self):
        match self.figs_row['plot_type']:
            case PlotType.LINE:
                self.line()
            case PlotType.SCATTER:
                self.scatter()
            case PlotType.ERRORBAR:
                NotImplementedError
            case PlotType.STEP:
                NotImplementedError
            case PlotType.FILLBETWEEN:
                NotImplementedError

            case PlotType.MESH2D:
                self.mesh2d()
            case PlotType.IMAGE2D:
                NotImplementedError
            case PlotType.CONTOUR2D:
                self.contour2d()
            case PlotType.CONTOURF2D:
                NotImplementedError
            case PlotType.MESHCONTOUR2D:
                NotImplementedError
            case PlotType.TRIMESH2D:
                NotImplementedError
            case PlotType.TRICONTOUR2D:
                NotImplementedError

            case PlotType.QUIVER2D:
                self.quiver2d()
            case PlotType.STREAM2D:
                self.stream2d()

            case PlotType.HIST1D:
                self.hist1d()
            case PlotType.HIST2D:
                self.hist2d()
            case PlotType.HEXBIN:
                NotImplementedError
            case PlotType.BOX:
                self.box()
            case PlotType.VIOLIN:
                self.violin()

            case PlotType.SURFACE3D:
                self.surface3d()
            case PlotType.WIRE3D:
                self.wire3d()
            case PlotType.SCATTER3D:
                self.scatter3d()
            case PlotType.LINE3D:
                NotImplementedError
            case PlotType.TRISURF3D:
                NotImplementedError
            case PlotType.CONTOUR3D:
                self.contour3d()
            case _:
                NotImplementedError

    def set_label_dict(self):
        self.label_dict = {}
        plot_type: str = self.figs_row['plot_type']
        
        if self.figs_row['legend_label'] is not None: self.label_dict['label'] = self.figs_row['legend_label']
        if plot_type=='line':
            color = self.figs_row['color']
            if color is not None:
                self.label_dict['color'] = color

    def set_common_axis_props(self):
        keys_and_functions: dict = {
            'xlabel': self.axis.set_xlabel,
            'xlim': self.axis.set_xlim,
            'xticks': self.axis.set_xticks,
            'xtick_labels': self.axis.set_xticklabels,
            'ylabel': self.axis.set_ylabel,
            'ylim': self.axis.set_ylim,
            'yticks': self.axis.set_yticks,
            'ytick_labels': self.axis.set_yticklabels,
            'title': self.axis.set_title,
        }

        if self.figs_row['plot_type'] in [
            PlotType.SURFACE3D,
            PlotType.WIRE3D   ,
            PlotType.SCATTER3D,
            PlotType.LINE3D   ,
            PlotType.TRISURF3D,
            PlotType.CONTOUR3D,
        ]:
            keys_and_functions.update({
                'zlabel': self.axis.set_zlabel,
                'zlim': self.axis.set_zlim,
                'zticks': self.axis.set_zticks,
                'ztick_labels': self.axis.set_zticklabels,
            })

            elev = self.figs_row['z_inc']
            azim = self.figs_row['z_azim']
            if elev is not None and azim is not None: 
                self.axis.view_init(elev=elev, azim=azim)
            
        for key, func in keys_and_functions.items():
            if self.figs_row[key] is not None: func(self.figs_row[key])

        legend_label = self.figs_row['legend_label']
        if legend_label is not None: self.axis.legend()

        xgrid: bool = self.figs_row['xgrid']
        ygrid: bool = self.figs_row['xgrid']
        if xgrid is not None: self.axis.grid(xgrid, axis='x')
        if ygrid is not None: self.axis.grid(ygrid, axis='y')
    
    def add_colorbar(self, mappable, label: str | None = None, **kwargs):
        fig: Figure = self.figs_row['figure']
        cbar = fig.colorbar(mappable, ax=self.axis, **kwargs)  # keeps it with this axes
        if label:
            cbar.set_label(label)
        return cbar

    def line(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_col: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        line_data: np.ndarray = self.dset_row['data'][data_col].to_numpy()

        self.set_label_dict()
        self.axis.plot(axis_data, line_data, **self.label_dict)
        self.set_common_axis_props()

    def scatter(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        scatter_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()

        additional_scatter_data = {}
        match scatter_data.ndim:
            case 1:
                np.expand_dims(scatter_data, axis=1)
            case 2:
                match scatter_data.shape[1]:
                    case 2:
                        additional_scatter_data['s'] = scatter_data[:, 1]
                    case 3:
                        additional_scatter_data['c'] = scatter_data[:, 2]
                        additional_scatter_data['cmap'] = 'viridis'
                    case _:
                        NotImplementedError

        self.set_label_dict()
        sc = self.axis.scatter(axis_data, scatter_data[:, 0], **additional_scatter_data, **self.label_dict)
        # If we used a colormap (i.e., 'c' in kwargs), add a colorbar:
        if 'c' in additional_scatter_data and 'cmap' in additional_scatter_data:
            self.add_colorbar(sc, label=self.figs_row.get('zlabel') or 'value')
        self.set_common_axis_props()

    def mesh2d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        mesh_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        Z = mesh_data.reshape(len(unique_x), len(unique_y))

        pcm = self.axis.pcolormesh(X, Y, Z, cmap='viridis', shading='gouraud')
        self.add_colorbar(pcm, label=self.figs_row.get('zlabel') or 'Z')
        self.set_common_axis_props()
    
    def contour2d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        contour_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        Z = contour_data.reshape(len(unique_x), len(unique_y))

        cf = self.axis.contourf(X, Y, Z, levels=20, cmap='viridis')
        self.add_colorbar(cf, label=self.figs_row.get('zlabel') or 'Z')
        self.set_common_axis_props()

    def quiver2d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        quiver_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        U = quiver_data[:, 0].reshape(len(unique_x), len(unique_y))
        V = quiver_data[:, 1].reshape(len(unique_x), len(unique_y))

        self.axis.quiver(X, Y, U, V)
        self.set_common_axis_props()

    def stream2d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        stream_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()

        # (Optional but recommended) make sure rows are ordered by x then y
        order = np.lexsort((axis_data[:, 1], axis_data[:, 0]))
        axis_data = axis_data[order]
        stream_data = stream_data[order]

        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])

        # Reshape to (Ny, Nx) => (len(unique_y), len(unique_x))
        U = stream_data[:, 0].reshape(len(unique_x), len(unique_y)).T
        V = stream_data[:, 1].reshape(len(unique_x), len(unique_y)).T

        # Pass 1-D axes directly
        self.axis.streamplot(unique_x, unique_y, U, V)
        self.set_common_axis_props()

    def hist1d(self):
        data_col: str = self.figs_row['dset_data_cols']
        hist_data: np.ndarray = self.dset_row['data'][data_col].to_numpy()

        self.set_label_dict()
        self.axis.hist(hist_data, bins=30, edgecolor='black', **self.label_dict)
        self.set_common_axis_props()

    def hist2d(self):
        data_col: str = self.figs_row['dset_data_cols']
        hist_data: np.ndarray = self.dset_row['data'][data_col].to_numpy()

        self.set_label_dict()
        H, xedges, yedges, im = self.axis.hist2d(hist_data[:, 0], hist_data[:, 1], bins=30, **self.label_dict)
        self.add_colorbar(im, label='count')
        self.set_common_axis_props()

    def box(self):
        data_cols: str = self.figs_row['dset_data_cols']
        box_df: pd.DataFrame = self.dset_row['data'][data_cols]
        box_data: np.ndarray = box_df.to_numpy()
        box_labels: list[str] = box_df.columns.to_list()

        self.set_label_dict()
        self.axis.boxplot(box_data)
        self.set_common_axis_props()

    def violin(self):
        data_cols: str = self.figs_row['dset_data_cols']
        violin_df: pd.DataFrame = self.dset_row['data'][data_cols]
        violin_data: np.ndarray = violin_df.to_numpy()
        violin_labels: list[str] = violin_df.columns.to_list()

        self.set_label_dict()
        self.axis.violinplot(violin_data)
        self.set_common_axis_props()

    def surface3d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        mesh_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        Z = mesh_data.reshape(len(unique_x), len(unique_y))

        surf = self.axis.plot_surface(X, Y, Z, cmap='viridis')
        self.add_colorbar(surf, label=self.figs_row.get('zlabel') or 'Z')
        self.set_common_axis_props()

    def wire3d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        mesh_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        Z = mesh_data.reshape(len(unique_x), len(unique_y))

        self.axis.plot_wireframe(X, Y, Z, cmap='viridis')
        self.set_common_axis_props()

    def scatter3d(self):
        axis_col: str = self.figs_row['dset_axis_cols']         # This will be x, y, z of the scatter point. 
        data_cols: str = self.figs_row['dset_data_cols']        # in this case data_cols will be color and size if included. 
        scatter_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        color_size_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()

        x = scatter_data[:, 0]
        y = scatter_data[:, 1]
        z = scatter_data[:, 2]

        color_size_data = {}
        if color_size_data is not None:
            match color_size_data.shape[1]:
                case 1:
                    color_size_data['s'] = color_size_data[:, 0]
                case 2:
                    color_size_data['color'] = color_size_data[:, 1]
                case _:
                    NotImplementedError

        self.axis.scatter(x, y, z, **color_size_data)
        self.set_common_axis_props()

    def contour3d(self):
        axis_col: str = self.figs_row['dset_axis_cols']
        data_cols: str = self.figs_row['dset_data_cols']
        axis_data: np.ndarray = self.dset_row['data'][axis_col].to_numpy()
        mesh_data: np.ndarray = self.dset_row['data'][data_cols].to_numpy()
        
        unique_x = np.unique(axis_data[:, 0])
        unique_y = np.unique(axis_data[:, 1])
        X, Y = np.meshgrid(unique_x, unique_y, indexing='ij')
        Z = mesh_data.reshape(len(unique_x), len(unique_y))

        self.axis.contour(X, Y, Z, zdir='z', offset=-1.0, cmap='viridis')
        self.set_common_axis_props()

class PlotType:
    # --- 2D basic ---
    LINE        = "line"
    SCATTER     = "scatter"
    ERRORBAR    = "errorbar"
    STEP        = "step"
    FILLBETWEEN = "fillbetween"

    # --- 2D field / image ---
    MESH2D          = "mesh2d"          # pcolormesh
    IMAGE2D         = "image2d"         # imshow
    CONTOUR2D       = "contour2d"       # contour lines
    CONTOURF2D      = "contourf2d"      # filled contour
    MESHCONTOUR2D   = "meshcontour2d"   # mesh + contour overlay
    TRIMESH2D       = "trimesh2d"       # tripcolor
    TRICONTOUR2D    = "tricontour2d"    # tricontour

    # --- 2D vector fields ---
    QUIVER2D    = "quiver2d"
    STREAM2D    = "stream2d"

    # --- distributions ---
    HIST1D      = "hist1d"
    HIST2D      = "hist2d"
    HEXBIN      = "hexbin"
    BOX         = "box"
    VIOLIN      = "violin"

    # --- 3D ---
    SURFACE3D   = "surface3d"
    WIRE3D      = "wire3d"
    SCATTER3D   = "scatter3d"
    LINE3D      = "line3d"
    TRISURF3D   = "trisurf3d"
    CONTOUR3D   = "contour3d"

    @classmethod
    def plot_requirements(cls):
        return {
            PlotType.LINE:        {"x": 1, "y": (1, None)},   # one x, one or many y
            PlotType.SCATTER:     {"x": 1, "y": (1, None)},
            PlotType.ERRORBAR:    {"x": 1, "y": 1, "yerr": (1, None)},
            PlotType.STEP:        {"x": 1, "y": (1, None)},
            PlotType.FILLBETWEEN: {"x": 1, "y": 2},           # y1 and y2 for shading

            PlotType.MESH2D:      {"x": 1, "y": 1, "z": 1},
            PlotType.IMAGE2D:     {"z": 1},
            PlotType.CONTOUR2D:   {"x": 1, "y": 1, "z": 1},
            PlotType.CONTOURF2D:  {"x": 1, "y": 1, "z": 1},
            PlotType.MESHCONTOUR2D: {"x": 1, "y": 1, "z": 1},
            PlotType.TRIMESH2D:   {"x": 1, "y": 1, "z": 1},
            PlotType.TRICONTOUR2D:{"x": 1, "y": 1, "z": 1},

            PlotType.QUIVER2D:    {"x": 1, "y": 1, "u": 1, "v": 1},
            PlotType.STREAM2D:    {"x": 1, "y": 1, "u": 1, "v": 1},

            PlotType.HIST1D:      {"x": 1},
            PlotType.HIST2D:      {"x": 1, "y": 1},
            PlotType.HEXBIN:      {"x": 1, "y": 1},
            PlotType.BOX:         {"y": (1, None)},
            PlotType.VIOLIN:      {"y": (1, None)},

            PlotType.SURFACE3D:   {"x": 1, "y": 1, "z": 1},
            PlotType.WIRE3D:      {"x": 1, "y": 1, "z": 1},
            PlotType.SCATTER3D:   {"x": 1, "y": 1, "z": 1},
            PlotType.LINE3D:      {"x": 1, "y": 1, "z": 1},
            PlotType.TRISURF3D:   {"x": 1, "y": 1, "z": 1},
            PlotType.CONTOUR3D:   {"x": 1, "y": 1, "z": 1},
        }

class PlotBase:
    def __init__(self, **kwargs):
        self.dsets_df: pd.DataFrame = pd.DataFrame({
            'name': pd.Series(dtype='string'),
            'data': pd.Series(dtype='object')
        })
        self.theme: str = 'default'
        self.figs_df: pd.DataFrame = pd.DataFrame({
            'fig_name': pd.Series(dtype='string'),
            'figure': pd.Series(dtype='object'),
            'subplot_nrow': pd.Series(dtype='Int64'),
            'subplot_ncol': pd.Series(dtype='Int64'),
            'subplot_idx': pd.Series(dtype='Int64'),
            'plot_type': pd.Series(dtype='string'),
            'axis': pd.Series(dtype='object'),
            'xlabel': pd.Series(dtype='string'),
            'xlim': pd.Series(dtype='object'),
            'xticks': pd.Series(dtype='object'),
            'xtick_labels': pd.Series(dtype='object'),
            'ylabel': pd.Series(dtype='string'),
            'ylim': pd.Series(dtype='object'),
            'yticks': pd.Series(dtype='object'),
            'ytick_labels': pd.Series(dtype='object'),
            'zlabel': pd.Series(dtype='string'),
            'zlim': pd.Series(dtype='object'),
            'zticks': pd.Series(dtype='object'),
            'ztick_labels': pd.Series(dtype='object'),
            'title': pd.Series(dtype='string'),
            'dset_name': pd.Series(dtype='string'),
            'dset_axis_cols': pd.Series(dtype='object'),
            'dset_data_cols': pd.Series(dtype='object'),
            'color': pd.Series(dtype='string'),
            'xgrid': pd.Series(dtype='bool'),
            'ygrid': pd.Series(dtype='bool'),
            'legend_label': pd.Series(dtype='string'),
        })

        self.inputdict: dict = None

        for key, value in kwargs.items():
            setattr(self, key, value)

    def save_data(self):
        os.makedirs('./plots', exist_ok=True)
        for idx, row in self.dsets_df.iterrows():
            name = row["name"]
            data = row['data']
            with h5py.File(f'./plots/{name}.h5', 'w') as f:
                f.create_dataset(f'{name}', shape=data.shape, data=data)

    def render_figures(self):
        set_common_rcparams(self.theme)
        
        # Set the figures and axis objects. 
        self.figs_cache = {}
        self.axis_cache = {}
        for idx, row in self.figs_df.iterrows():
            fig_name = row["fig_name"]
            color = row['color']
            plot_type = row['plot_type']
            plots_that_are_3d = [
                PlotType.SURFACE3D,
                PlotType.WIRE3D,
                PlotType.SCATTER3D,
                PlotType.LINE3D,
                PlotType.TRISURF3D,
                PlotType.CONTOUR3D,
            ]
            kwargs_dict = {}
            if plot_type in plots_that_are_3d: kwargs_dict['projection'] = '3d'

            # ---- FIGURE CACHE ----
            if fig_name not in self.figs_cache:
                self.figs_cache[fig_name] = plt.figure()
            figure = self.figs_cache[fig_name]
            self.figs_df.at[idx, "figure"] = figure

            # ---- AXIS CACHE ----
            axis_key = (fig_name, row["subplot_nrow"], row["subplot_ncol"], row["subplot_idx"])
            if axis_key not in self.axis_cache:
                self.axis_cache[axis_key] = figure.add_subplot(
                    row["subplot_nrow"], row["subplot_ncol"], row["subplot_idx"], **kwargs_dict
                )
            axis = self.axis_cache[axis_key]
            self.figs_df.at[idx, "axis"] = axis

        # plot.
        for idx, fig_row in self.figs_df.iterrows():
            axis = fig_row['axis']
            dset_name = fig_row['dset_name']
            dset_row = self.dsets_df.loc[self.dsets_df['name'] == dset_name].iloc[0]
            axisplot = AxisPlot(
                axis=axis,
                dset_row=dset_row,
                figs_row=fig_row,
                dsets_df=self.dsets_df,
                figs_df=self.figs_df,
            )
            axisplot.plot()

    def save_figures(self, **kwargs):
        os.makedirs('./plots', exist_ok=True)
        # Render figures. 
        self.render_figures()

        # Save them. 
        for key, group in self.figs_df[['fig_name', 'figure']].groupby(by=['fig_name']):
            figname = key[0]
            fig: Figure = group['figure'].iloc[0]
            fig.savefig(f'./plots/{figname}.png')

#endregion