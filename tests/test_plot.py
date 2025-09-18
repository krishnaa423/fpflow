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
from fpflow.plots.plot import PlotBase, PlotType
#endregion

#region variables
#endregion

#region functions
def main():
    tp = TestPlot()
    tp.save_figures()

#endregion

#region classes
class TestPlot(PlotBase):
    def __init__(self):
        super().__init__()

        # -----------------------------
        # Build datasets to cover cases
        # -----------------------------
        # 1) Simple band-like lines
        dset_bands = pd.DataFrame({
            'x': np.linspace(0, 1, 200),
            'y1': np.sin(2*np.pi*np.linspace(0,1,200)) + 0.2*np.random.randn(200),
            'y2': np.cos(2*np.pi*np.linspace(0,1,200)) + 0.2*np.random.randn(200),
            'y3': np.linspace(0, 1, 200),
            'y4': np.linspace(1, 0, 200),
        })

        # 2) Scatter with marker sizes (2 cols => size used)
        n_sc = 250
        x_sc = np.linspace(-2, 2, n_sc)
        y_sc = np.sinc(x_sc) + 0.05*np.random.randn(n_sc)
        s_sc = 20 + 80*(np.abs(x_sc)/x_sc.max())  # size
        dset_scatter_size = pd.DataFrame({'x': x_sc, 'y': y_sc, 's': s_sc})

        # 3) Scatter with colors (3 cols => color used)
        n_sc2 = 300
        x2 = np.linspace(-3, 3, n_sc2)
        y2 = np.sin(x2) + 0.1*np.random.randn(n_sc2)
        c2 = x2  # color by x
        dset_scatter_color = pd.DataFrame({'x': x2, 'y': y2, 's': 30*np.ones(n_sc2), 'c': c2})

        # 4) Regular grid for field plots (mesh/contour/quiver/stream + 3D)
        nx, ny = 40, 25
        gx = np.linspace(-2*np.pi, 2*np.pi, nx)
        gy = np.linspace(-1.5*np.pi, 1.5*np.pi, ny)
        GX, GY = np.meshgrid(gx, gy, indexing='xy')
        Z = np.sin(GX)*np.cos(GY)

        # Vector field (a simple rotational field)
        U = -np.gradient(Z, axis=1)  # d/dy (approx)
        V =  np.gradient(Z, axis=0)  # d/dx (approx)

        # Flatten to match AxisPlot’s reshape(len(unique_x), len(unique_y)) assumption
        XY_flat = np.column_stack([GX.ravel(), GY.ravel()])
        Z_flat  = Z.ravel()
        U_flat  = U.ravel()
        V_flat  = V.ravel()

        dset_grid_scalar = pd.DataFrame({'x': XY_flat[:,0], 'y': XY_flat[:,1], 'z': Z_flat})
        dset_grid_vector = pd.DataFrame({'x': XY_flat[:,0], 'y': XY_flat[:,1], 'u': U_flat, 'v': V_flat})

        # 5) Distributions
        rng = np.random.default_rng(7)
        dset_hist1d = pd.DataFrame({'samples': rng.normal(loc=0.0, scale=1.0, size=2000)})

        # Correlated 2D normal for hist2d
        mean = np.array([0.0, 0.0])
        cov = np.array([[1.0, 0.75],
                        [0.75, 1.5]])
        xy = rng.multivariate_normal(mean, cov, size=2500)
        dset_hist2d = pd.DataFrame({'x': xy[:,0], 'y': xy[:,1]})

        # 6) Box/Violin data
        dset_boxviolin = pd.DataFrame({
            'A': rng.normal(0.0, 1.0, 400),
            'B': rng.normal(1.0, 0.6, 400),
            'C': rng.normal(-0.5, 1.2, 400),
            'D': rng.uniform(-2.0, 2.0, 400),
        })

        # 7) 3D surface uses same grid scalar
        #    (x,y,z columns already present in dset_grid_scalar)

        dsets_append: pd.DataFrame = pd.DataFrame([
            {'name': 'dset_bands',          'data': dset_bands},
            {'name': 'dset_scatter_size',   'data': dset_scatter_size},
            {'name': 'dset_scatter_color',  'data': dset_scatter_color},
            {'name': 'dset_grid_scalar',    'data': dset_grid_scalar},
            {'name': 'dset_grid_vector',    'data': dset_grid_vector},
            {'name': 'dset_hist1d',         'data': dset_hist1d},
            {'name': 'dset_hist2d',         'data': dset_hist2d},
            {'name': 'dset_boxviolin',      'data': dset_boxviolin},
        ])

        # -----------------------------------
        # Figures to hit every implemented API
        # -----------------------------------
        figs_append: pd.DataFrame = pd.DataFrame([
            # ---- LINE ----
            {
                'fig_name': 'fig_line',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 2, 'subplot_idx': 1,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'k-point', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'Energy (a.u.)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Line Plot',
                'dset_name': 'dset_bands',
                'dset_axis_cols': 'x',          # IMPORTANT: use string here
                'dset_data_cols': 'y1',
                'color': None, 'legend_label': 'DFT',
            },
            {
                'fig_name': 'fig_line',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 2, 'subplot_idx': 2,
                'plot_type': PlotType.LINE, 'axis': None,
                'xlabel': 'k-point', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'Energy (a.u.)', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Line Plot',
                'dset_name': 'dset_bands',
                'dset_axis_cols': 'x',          # IMPORTANT: use string here
                'dset_data_cols': 'y2',
                'color': None, 'legend_label': 'DFT',
            },
            # ---- SCATTER (size) ----
            {
                'fig_name': 'fig_scatter_size',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTER, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Scatter (size)',
                'dset_name': 'dset_scatter_size',
                'dset_axis_cols': 'x',
                'dset_data_cols': ['y', 's'],   # 2 cols => size is used
                'color': None, 'legend_label': 'points',
            },
            # ---- SCATTER (color) ----
            {
                'fig_name': 'fig_scatter_color',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SCATTER, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Scatter (color)',
                'dset_name': 'dset_scatter_color',
                'dset_axis_cols': 'x',
                'dset_data_cols': ['y', 's', 'c'],  # 3 cols => color is used (per your method)
                'color': None, 'legend_label': 'colored',
            },
            # ---- MESH2D ----
            {
                'fig_name': 'fig_mesh2d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.MESH2D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'pcolormesh',
                'dset_name': 'dset_grid_scalar',
                'dset_axis_cols': ['x', 'y'],   # list -> 2D axis array
                'dset_data_cols': 'z',
                'color': None, 'legend_label': None,
            },
            # ---- CONTOUR2D ----
            {
                'fig_name': 'fig_contour2d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.CONTOUR2D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Filled Contour',
                'dset_name': 'dset_grid_scalar',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': 'z',
                'color': None, 'legend_label': None,
            },
            # ---- QUIVER2D ----
            {
                'fig_name': 'fig_quiver2d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.QUIVER2D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Quiver',
                'dset_name': 'dset_grid_vector',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': ['u', 'v'],
                'color': None, 'legend_label': None,
            },
            # ---- STREAM2D ----
            {
                'fig_name': 'fig_stream2d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.STREAM2D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Streamplot',
                'dset_name': 'dset_grid_vector',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': ['u', 'v'],
                'color': None, 'legend_label': None,
            },
            # ---- HIST1D ----
            {
                'fig_name': 'fig_hist1d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.HIST1D, 'axis': None,
                'xlabel': 'value', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'count', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Histogram 1D',
                'dset_name': 'dset_hist1d',
                'dset_axis_cols': 'samples',   # not used, but harmless
                'dset_data_cols': 'samples',
                'color': None, 'legend_label': 'N(0,1)',
            },
            # ---- HIST2D ----
            {
                'fig_name': 'fig_hist2d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.HIST2D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Histogram 2D',
                'dset_name': 'dset_hist2d',
                'dset_axis_cols': ['x', 'y'],   # not used by hist2d(), but safe
                'dset_data_cols': ['x', 'y'],   # EXPECTED by your hist2d()
                'color': None, 'legend_label': None,
            },
            # ---- BOX ----
            {
                'fig_name': 'fig_box',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.BOX, 'axis': None,
                'xlabel': None, 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': None, 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Boxplot',
                'dset_name': 'dset_boxviolin',
                'dset_axis_cols': None,
                'dset_data_cols': ['A', 'B', 'C', 'D'],
                'color': None, 'legend_label': None,
            },
            # ---- VIOLIN ----
            {
                'fig_name': 'fig_violin',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.VIOLIN, 'axis': None,
                'xlabel': None, 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': None, 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': None, 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': None, 'z_azim': None,
                'title': 'Violin',
                'dset_name': 'dset_boxviolin',
                'dset_axis_cols': None,
                'dset_data_cols': ['A', 'B', 'C', 'D'],
                'color': None, 'legend_label': None,
            },
            # ---- SURFACE3D ----
            {
                'fig_name': 'fig_surface3d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.SURFACE3D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': 'z', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': 25, 'z_azim': -60,
                'title': 'Surface 3D',
                'dset_name': 'dset_grid_scalar',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': 'z',
                'color': None, 'legend_label': None,
            },
            # ---- WIRE3D ----
            {
                'fig_name': 'fig_wire3d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.WIRE3D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': 'z', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': 30, 'z_azim': 35,
                'title': 'Wireframe 3D',
                'dset_name': 'dset_grid_scalar',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': 'z',
                'color': None, 'legend_label': None,
            },
            # ---- CONTOUR3D ----
            {
                'fig_name': 'fig_contour3d',
                'figure': None, 'subplot_nrow': 1, 'subplot_ncol': 1, 'subplot_idx': 1,
                'plot_type': PlotType.CONTOUR3D, 'axis': None,
                'xlabel': 'x', 'xlim': None, 'xticks': None, 'xtick_labels': None,
                'ylabel': 'y', 'ylim': None, 'yticks': None, 'ytick_labels': None,
                'zlabel': 'z', 'zlim': None, 'zticks': None, 'ztick_labels': None,
                'z_inc': 35, 'z_azim': -30,
                'title': 'Contour 3D (offset in z)',
                'dset_name': 'dset_grid_scalar',
                'dset_axis_cols': ['x', 'y'],
                'dset_data_cols': 'z',
                'color': None, 'legend_label': None,
            },
        ])

        # Apply
        self.dsets_df = pd.concat([self.dsets_df, dsets_append], ignore_index=True)
        self.figs_df  = pd.concat([self.figs_df,  figs_append],  ignore_index=True)

#endregion

#region main
if __name__=='__main__':
    main()
#endregion