# plot_utils.py
# All the functions to do the plots for the manuscripts

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import nibabel as nib
import numpy as np
import pandas as pd
from neuromaps.datasets import fetch_fslr

from nilearn import surface, plotting, datasets, image

######### the standard surface maps ##################
def plot_surface_maps(
    left_gii,
    right_gii,
    out_file,
    colorbar_label="Effect size",
    current_threshold=None
):
    """
    Plot a 4-panel cortical surface map:
    left/right hemisphere × lateral/medial views.

    Parameters
    ----------
    left_gii : str
        Path to left hemisphere GIFTI file.
    right_gii : str
        Path to right hemisphere GIFTI file.
    out_file : str
        Output image file path.
    colorbar_label : str, optional
        Label displayed on the colorbar.
    """

    # ====================================
    # 1. Figure style
    # ====================================
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # ====================================
    # 2. Load surface data
    # ====================================
    texture_left = nib.load(left_gii).darrays[0].data
    texture_right = nib.load(right_gii).darrays[0].data

    # ====================================
    # 3. Load fsaverage surfaces
    # ====================================
    fsaverage = datasets.fetch_surf_fsaverage()

    # ====================================
    # 4. Define color scale
    # ====================================
    all_data = np.concatenate([texture_left, texture_right])

    absmax = np.nanmax(np.abs(all_data))
    vmin, vmax = -absmax, absmax

    cmap = "seismic"
    norm = TwoSlopeNorm(
        vmin=vmin,
        vcenter=0,
        vmax=vmax
    )

    # ====================================
    # 5. Create figure
    # ====================================
    fig = plt.figure(figsize=(12, 7), dpi=300)

    axes = [
        fig.add_subplot(2, 2, 1, projection="3d"),
        fig.add_subplot(2, 2, 2, projection="3d"),
        fig.add_subplot(2, 2, 3, projection="3d"),
        fig.add_subplot(2, 2, 4, projection="3d"),
    ]

    # ====================================
    # 6. Plot surfaces
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_left,
        texture_left,
        hemi="left",
        view="lateral",
        bg_map=fsaverage.sulc_left,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[0],
        threshold=current_threshold,
    )

    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture_right,
        hemi="right",
        view="lateral",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[1],
        threshold=current_threshold
    )

    plotting.plot_surf_stat_map(
        fsaverage.pial_left,
        texture_left,
        hemi="left",
        view="medial",
        bg_map=fsaverage.sulc_left,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[2],
        threshold=current_threshold
    )

    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture_right,
        hemi="right",
        view="medial",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[3],
        threshold=current_threshold
    )

    # ====================================
    # 7. Enlarge brain panels and reduce whitespace
    # ====================================
    
    axes[0].set_position([0.01, 0.58, 0.58, 0.58])   # LH lateral
    axes[1].set_position([0.31, 0.58, 0.58, 0.58])   # RH lateral
    
    axes[2].set_position([0.01, 0.2, 0.58, 0.58])   # LH medial
    axes[3].set_position([0.31, 0.2, 0.58, 0.58])   # RH medial
    
    for ax in axes:
        ax.set_facecolor((0, 0, 0, 0))
    
   # ====================================
    # 8. Horizontal colorbar
    # ====================================
    sm = ScalarMappable(
        norm=norm,
        cmap=cmap
    )
    sm.set_array([])
    
    cbar_ax = fig.add_axes(
        [0.25, 0.23, 0.40, 0.035]
    )
    
    cbar = fig.colorbar(
        sm,
        cax=cbar_ax,
        orientation="horizontal"
    )
    
    # label above colorbar
    cbar.set_label(
        colorbar_label,
        fontsize=24,
        fontweight="bold",
        labelpad=8
    )
    
    cbar.ax.xaxis.set_label_position("top")
    
    # ticks
    cbar.ax.tick_params(
        axis="x",
        which="both",
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True,
        labelsize=16,
        length=0
    )
    
    for tick in cbar.ax.get_xticklabels():
        tick.set_fontweight("bold")
    
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    
    # remove border
    cbar.outline.set_visible(False)

    # ====================================
    # 9. Save figure
    # ====================================
    os.makedirs(
        os.path.dirname(out_file),
        exist_ok=True
    )

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
        transparent=True
    )

    plt.show()

    return out_file



############# four surface views in one row ##################
def plot_surface_maps_onerow(
    left_gii,
    right_gii,
    out_file,
    colorbar_label="Effect size",
    current_threshold=None
):
    """Plot LH lateral, LH medial, RH lateral, RH medial in one row.

    Accepts the same arguments as ``plot_surface_maps`` and uses its
    surface, font, color scale, and export style, with a vertical colorbar
    at the far right. Returns the output image path.
    """
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    texture_left = nib.load(left_gii).darrays[0].data
    texture_right = nib.load(right_gii).darrays[0].data
    fsaverage = datasets.fetch_surf_fsaverage()

    all_data = np.concatenate([texture_left, texture_right])
    finite_data = all_data[np.isfinite(all_data)]
    if finite_data.size == 0:
        raise ValueError("Surface maps must contain at least one finite value.")
    # Thresholded maps can contain no surviving vertices.
    absmax = np.max(np.abs(finite_data))
    if absmax == 0:
        absmax = 1.0
    vmin, vmax = -absmax, absmax
    cmap = "seismic"
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(18, 4), dpi=300)
    views = [
        ("left", "lateral", fsaverage.pial_left, texture_left, fsaverage.sulc_left),
        ("left", "medial", fsaverage.pial_left, texture_left, fsaverage.sulc_left),
        ("right", "lateral", fsaverage.pial_right, texture_right, fsaverage.sulc_right),
        ("right", "medial", fsaverage.pial_right, texture_right, fsaverage.sulc_right),
    ]
    for i, (hemi, view, mesh, texture, sulc) in enumerate(views):
        ax = fig.add_subplot(1, 4, i + 1, projection="3d")
        plotting.plot_surf_stat_map(
            mesh,
            texture,
            hemi=hemi,
            view=view,
            bg_map=sulc,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            colorbar=False,
            axes=ax,
            threshold=current_threshold,
        )
        # Overlap the axes' empty margins to reduce visible brain gaps by ~2/3.
        ax.set_position([0.005 + i * 0.188, 0.02, 0.225, 0.96])
        ax.set_facecolor((0, 0, 0, 0))

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.809, 0.20, 0.012, 0.60])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
    cbar.set_label(colorbar_label, fontsize=24, fontweight="bold", labelpad=8)
    cbar.ax.yaxis.set_label_position("right")
    cbar.ax.tick_params(
        axis="y", which="both", left=False, right=True,
        labelleft=False, labelright=True, labelsize=16, length=0
    )
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight("bold")
    cbar.outline.set_visible(False)

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    fig.savefig(
        out_file, dpi=300, bbox_inches="tight", pad_inches=0, transparent=True
    )
    plt.show()
    return out_file


############# the single column surface maps ##################
def plot_surface_maps_singlecol(
    left_gii,
    right_gii,
    out_file,
    colorbar_label="Effect size",
    current_threshold=None
):
    """
    Plot a 2-panel cortical surface map:
    left medial view + right lateral view

    Parameters
    ----------
    left_gii : str
        Path to left hemisphere GIFTI file.
    right_gii : str
        Path to right hemisphere GIFTI file.
    out_file : str
        Output image file path.
    colorbar_label : str, optional
        Label displayed on the colorbar.
    """

    # ====================================
    # 1. Figure style
    # ====================================
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # ====================================
    # 2. Load surface data
    # ====================================
    texture_left = nib.load(left_gii).darrays[0].data
    texture_right = nib.load(right_gii).darrays[0].data

    # ====================================
    # 3. Load fsaverage surfaces
    # ====================================
    fsaverage = datasets.fetch_surf_fsaverage()

    # ====================================
    # 4. Define color scale
    # ====================================
    all_data = np.concatenate([texture_left, texture_right])

    absmax = np.nanmax(np.abs(all_data))
    vmin, vmax = -absmax, absmax

    cmap = "seismic"
    norm = TwoSlopeNorm(
        vmin=vmin,
        vcenter=0,
        vmax=vmax
    )

    # ====================================
    # 5. Create figure
    # ====================================
    fig = plt.figure(figsize=(4.5, 6), dpi=300)

    axes = [
        fig.add_subplot(2, 1, 1, projection="3d"),
        fig.add_subplot(2, 1, 2, projection="3d")
    ]

    # ====================================
    # 6. Plot surfaces
    # ====================================

    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture_right,
        hemi="right",
        view="lateral",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[0],
        threshold=current_threshold
    )
    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture_right,
        hemi="right",
        view="medial",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        axes=axes[1],
        threshold=current_threshold
    )

    # ====================================
    # 7. Enlarge brain panels and reduce whitespace
    # ====================================
    
    axes[0].set_position([0.00, 0.46, 0.92, 0.52])
    axes[1].set_position([0.00, 0.12, 0.92, 0.52])

    for ax in axes:
        ax.set_facecolor((0, 0, 0, 0))
    
    # ====================================
    # 8. Vertical colorbar
    # ====================================
    
    sm = ScalarMappable(
        norm=norm,
        cmap=cmap
    )
    sm.set_array([])
    
    cbar_ax = fig.add_axes(
        [0.85, 0.25, 0.04, 0.60]
    )
    
    cbar = fig.colorbar(
        sm,
        cax=cbar_ax,
        orientation="vertical"
    )
    
    cbar.set_label(
        colorbar_label,
        fontsize=18,
        fontweight="bold",
        rotation=90,
        labelpad=5
    )
    
    cbar.ax.yaxis.set_label_position("left")
    
    cbar.ax.tick_params(
        labelsize=14,
        length=0
    )
    
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight("bold")
    
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    
    cbar.outline.set_visible(False)

    # ====================================
    # 9. Save figure
    # ====================================
    os.makedirs(
        os.path.dirname(out_file),
        exist_ok=True
    )

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
        transparent=True
    )

    plt.show()

    return out_file


#################################################################
######### Right hemisphere: medial + lateral (one row) ##########
#################################################################

def plot_surface_maps_singlerow(
    right_gii,
    out_file,
    current_threshold=None,
    current_cmap = "seismic",
):

    # ====================================
    # Figure style
    # ====================================
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # ====================================
    # Load data
    # ====================================
    texture = nib.load(right_gii).darrays[0].data

    # ====================================
    # fsaverage
    # ====================================
    fsaverage = datasets.fetch_surf_fsaverage()
    #surfaces = fetch_fslr(density="32k")

    # ====================================
    # Color scale
    # ====================================
    absmax = np.nanmax(np.abs(texture))
    vmin = -absmax
    vmax = absmax

    cmap = current_cmap

    # ====================================
    # Figure
    # ====================================
    fig = plt.figure(figsize=(8, 3.8), dpi=300)

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    # ====================================
    # Right medial
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        #surfaces["inflated"][1],
        texture,
        hemi="right",
        view="medial",
        bg_map=fsaverage.sulc_right,
        #bg_map=surfaces["sulc"][1],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax1,
    )

    # ====================================
    # Right lateral
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        #surfaces["inflated"][1],
        texture,
        hemi="right",
        view="lateral",
        bg_map=fsaverage.sulc_right,
        #bg_map=surfaces["sulc"][1],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax2,
    )

    # ====================================
    # Adjust layout
    # ====================================
    ax1.set_position([-0.01, -0.02, 0.58, 1.04])
    ax2.set_position([0.43, -0.02, 0.58, 1.04])

    for ax in [ax1, ax2]:
        ax.set_facecolor((0, 0, 0, 0))

    # ====================================
    # Save
    # ====================================
    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
        transparent=True,
    )
    plt.show()
    return out_file


######### Left hemisphere parcels: medial + lateral (one row) ##########
def plot_surface_pacel_lh(
    parcel_vector,
    out_file,
    current_threshold=None,
    label_gii=(
        "/Users/Chenxiao/Documents/MATLAB/DPABI_V9.0_250415/"
        "DPABISurf/SurfTemplates/"
        "fsaverage5_lh_Schaefer2018_400Parcels_7Networks_order.label.gii"
    ),
    current_cmap = "seismic",
):
    """
    Map 200 left-hemisphere Schaefer parcel values to fsaverage5 vertices
    and plot medial and lateral views in one row.

    Parameters
    ----------
    parcel_vector : array-like
        A vector containing exactly 200 values. Position 1 maps to label 1,
        position 2 maps to label 2, and so on through label 200.
    out_file : str
        Output image file path.
    current_threshold : float, optional
        Absolute plotting threshold passed to ``plot_surf_stat_map``.
    label_gii : str, optional
        Path to the left-hemisphere Schaefer-400 label GIFTI. The default is
        the DPABISurf fsaverage5 template specified for this project.

    Returns
    -------
    str
        The output image path.
    """

    # ====================================
    # Figure style
    # ====================================
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # ====================================
    # Map parcel values to vertices
    # ====================================
    parcel_values = np.asarray(parcel_vector, dtype=float).reshape(-1)
    if parcel_values.size != 200:
        raise ValueError(
            "parcel_vector must contain exactly 200 values for the "
            "left-hemisphere parcels."
        )
    if not np.isfinite(parcel_values).any():
        raise ValueError("parcel_vector must contain at least one finite value.")

    label_img = nib.load(label_gii)
    parcel_labels = np.asarray(label_img.darrays[0].data).squeeze()
    if parcel_labels.ndim != 1:
        raise ValueError("The parcel label GIFTI must contain one-dimensional data.")

    expected_labels = np.arange(1, 201)
    observed_labels = np.unique(parcel_labels[parcel_labels > 0])
    if not np.array_equal(observed_labels, expected_labels):
        raise ValueError(
            "The parcel label GIFTI must contain cortical labels 1 through 200."
        )

    texture = np.zeros(parcel_labels.shape, dtype=float)
    cortex_mask = parcel_labels > 0
    texture[cortex_mask] = parcel_values[parcel_labels[cortex_mask].astype(int) - 1]

    # ====================================
    # fsaverage
    # ====================================
    fsaverage = datasets.fetch_surf_fsaverage(mesh="fsaverage5")
    n_surface_vertices = surface.load_surf_mesh(
        fsaverage.pial_left
    ).coordinates.shape[0]
    if texture.size != n_surface_vertices:
        raise ValueError(
            "The label GIFTI and fsaverage5 left surface have different "
            f"vertex counts: {texture.size} versus {n_surface_vertices}."
        )

    # ====================================
    # Color scale
    # ====================================
    finite_texture = texture[np.isfinite(texture)]
    absmax = np.max(np.abs(finite_texture))
    if absmax == 0:
        absmax = 1.0

    vmin = -absmax
    vmax = absmax
    cmap = current_cmap

    # ====================================
    # Figure
    # ====================================
    fig = plt.figure(figsize=(8, 3.8), dpi=300)

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    # ====================================
    # Left medial
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_left,
        texture,
        hemi="left",
        view="medial",
        bg_map=fsaverage.sulc_left,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax1,
    )

    # ====================================
    # Left lateral
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_left,
        texture,
        hemi="left",
        view="lateral",
        bg_map=fsaverage.sulc_left,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax2,
    )

    # ====================================
    # Adjust layout
    # ====================================
    ax1.set_position([-0.01, -0.02, 0.58, 1.04])
    ax2.set_position([0.43, -0.02, 0.58, 1.04])

    for ax in [ax1, ax2]:
        ax.set_facecolor((0, 0, 0, 0))

    # ====================================
    # Save
    # ====================================
    output_dir = os.path.dirname(os.path.abspath(out_file))
    os.makedirs(output_dir, exist_ok=True)

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
        transparent=True,
    )

    plt.show()

    return out_file


######### Right hemisphere parcels: medial + lateral (one row) ##########
def plot_surface_pacel_rh(
    parcel_vector,
    out_file,
    current_threshold=None,
    label_gii=(
        "/Users/Chenxiao/Documents/MATLAB/DPABI_V9.0_250415/"
        "DPABISurf/SurfTemplates/"
        "fsaverage5_rh_Schaefer2018_400Parcels_7Networks_order.label.gii"
    ),
    current_cmap="Oranges",
    current_vmin=None,
    current_vmax=None,
):
    """
    Map 200 right-hemisphere Schaefer parcel values to fsaverage5 vertices
    and plot medial and lateral views in one row.

    Parameters
    ----------
    parcel_vector : array-like
        A vector containing exactly 200 values. Position 1 maps to label 1,
        position 2 maps to label 2, and so on through label 200.
    out_file : str
        Output image file path.
    current_threshold : float, optional
        Absolute plotting threshold passed to ``plot_surf_stat_map``.
    label_gii : str, optional
        Path to the right-hemisphere Schaefer-400 label GIFTI. The default is
        the DPABISurf fsaverage5 template specified for this project.
    current_cmap : str or matplotlib.colors.Colormap, optional
        Colormap used for the parcel values. Defaults to ``"Oranges"``.

    current_vmin, current_vmax : float, optional
        Explicit color limits. Omitted limits retain the symmetric defaults.
        For nonnegative expression, use zero and the finite parcel maximum.

    Returns
    -------
    str
        The output image path.
    """

    # ====================================
    # Figure style
    # ====================================
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # ====================================
    # Map parcel values to vertices
    # ====================================
    parcel_values = np.asarray(parcel_vector, dtype=float).reshape(-1)
    if parcel_values.size != 200:
        raise ValueError(
            "parcel_vector must contain exactly 200 values for the "
            "right-hemisphere parcels."
        )
    if not np.isfinite(parcel_values).any():
        raise ValueError("parcel_vector must contain at least one finite value.")

    label_img = nib.load(label_gii)
    parcel_labels = np.asarray(label_img.darrays[0].data).squeeze()
    if parcel_labels.ndim != 1:
        raise ValueError("The parcel label GIFTI must contain one-dimensional data.")

    expected_labels = np.arange(1, 201)
    observed_labels = np.unique(parcel_labels[parcel_labels > 0])
    if not np.array_equal(observed_labels, expected_labels):
        raise ValueError(
            "The parcel label GIFTI must contain cortical labels 1 through 200."
        )

    texture = np.zeros(parcel_labels.shape, dtype=float)
    cortex_mask = parcel_labels > 0
    texture[cortex_mask] = parcel_values[parcel_labels[cortex_mask].astype(int) - 1]

    # ====================================
    # fsaverage
    # ====================================
    fsaverage = datasets.fetch_surf_fsaverage(mesh="fsaverage5")
    n_surface_vertices = surface.load_surf_mesh(
        fsaverage.pial_right
    ).coordinates.shape[0]
    if texture.size != n_surface_vertices:
        raise ValueError(
            "The label GIFTI and fsaverage5 right surface have different "
            f"vertex counts: {texture.size} versus {n_surface_vertices}."
        )

    # ====================================
    # Color scale
    # ====================================
    finite_texture = texture[np.isfinite(texture)]
    absmax = np.max(np.abs(finite_texture))
    if absmax == 0:
        absmax = 1.0

    vmin = -absmax if current_vmin is None else float(current_vmin)
    vmax = absmax if current_vmax is None else float(current_vmax)
    if not np.isfinite([vmin, vmax]).all() or vmin >= vmax:
        raise ValueError("Color limits must be finite and satisfy vmin < vmax.")
    cmap = current_cmap

    # ====================================
    # Figure
    # ====================================
    fig = plt.figure(figsize=(8, 3.8), dpi=300)

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    # ====================================
    # Right medial
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture,
        hemi="right",
        view="medial",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax1,
    )

    # ====================================
    # Right lateral
    # ====================================
    plotting.plot_surf_stat_map(
        fsaverage.pial_right,
        texture,
        hemi="right",
        view="lateral",
        bg_map=fsaverage.sulc_right,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        colorbar=False,
        threshold=current_threshold,
        axes=ax2,
    )

    # ====================================
    # Adjust layout
    # ====================================
    ax1.set_position([-0.01, -0.02, 0.58, 1.04])
    ax2.set_position([0.43, -0.02, 0.58, 1.04])

    for ax in [ax1, ax2]:
        ax.set_facecolor((0, 0, 0, 0))

    # ====================================
    # Save
    # ====================================
    output_dir = os.path.dirname(os.path.abspath(out_file))
    os.makedirs(output_dir, exist_ok=True)

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0,
        transparent=True,
    )

    plt.show()

    return out_file


######################################################   
##########          radar plot           #############
######################################################
def plot_radar(
        csv_file,
        output_file="radar_plot.png",
        col1="Group1",
        col2="Group2",
        label1="Group1",
        label2="Group2",
        color1="#1565C0",
        color2="#C62828",
        dpi=600,
        font_size=18,
        figsize=(6, 6),
        legend_fontsize=14,
        tick_fontsize=14):
    """Plot two groups as a radar chart.

    font_size controls network labels; tick_fontsize controls radial tick labels.
    Font sizes (including legend_fontsize) are in points. figsize is (width,
    height) in inches, before the tight bounding box used when saving.
    """

    # ---------------------
    # Read data
    # ---------------------
    df = pd.read_csv(csv_file)

    labels = df["Network"].tolist()

    values1 = df[col1].values
    values2 = df[col2].values

    # Close polygon
    values1 = np.append(values1, values1[0])
    values2 = np.append(values2, values2[0])

    angles = np.linspace(
        0,
        2 * np.pi,
        len(labels),
        endpoint=False
    )

    angles = np.append(angles, angles[0])

    # ---------------------
    # Figure style
    # ---------------------
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 14,
        "axes.linewidth": 1.5
    })

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111, polar=True)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Grid
    ax.grid(
        linestyle="--",
        linewidth=0.8,
        alpha=0.5
    )

    # Labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=font_size) #fontweight='bold'
    ax.tick_params(axis='x',pad=24)

    yticks = ax.get_yticks()
    ax.set_yticks(yticks)
    ax.set_yticklabels(
    [f"{y:.1f}" for y in yticks[:-1]] + [''], fontsize=tick_fontsize
    )

    # Lines
    ax.plot(
        angles,
        values1,
        color=color1,
        linewidth=2.5,
        label=label1
    )

    ax.plot(
        angles,
        values2,
        color=color2,
        linewidth=2.5,
        label=label2
    )

    # Fill
    ax.fill(
        angles,
        values1,
        color=color1,
        alpha=0.15
    )

    ax.fill(
        angles,
        values2,
        color=color2,
        alpha=0.15
    )

    ax.spines["polar"].set_linewidth(1)

    ax.legend(
        frameon=False,
        fontsize=legend_fontsize,
        bbox_to_anchor=(1.2, 1.1)
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        transparent=True
    )

    plt.show()

    print(f"Saved: {output_file}")
