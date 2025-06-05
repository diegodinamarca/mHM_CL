# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to visualize LAI preprocessing for quality control.
Generates maps of original, filled, and reprojected LAI rasters, and maps from the NetCDF file.
"""

import os
import json
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

out_folder = config["out_folder"]
lai_folder = config["lai_folder"]

output_folder = os.path.join(out_folder, "Preprocess LAI QC")
os.makedirs(output_folder, exist_ok=True)

lai_temp_folder = os.path.join(lai_folder, "temp")

# ==== File Paths ====

clipped_lai_fp = os.path.join(lai_temp_folder, "clipped_lai.tif")
clipped_filled_lai_fp = os.path.join(lai_temp_folder, "clipped_filled_lai.tif")
reprojected_lai_fp = os.path.join(lai_temp_folder, "reprojected_lai.tif")
lai_nc_fp = os.path.join(lai_folder, "lai.nc")

# ==== Helper Functions ====

def load_raster_band(fp, band=1):
    with rasterio.open(fp) as src:
        arr = src.read(band)
        arr = np.where(arr == src.nodata, np.nan, arr)
    return arr

def plot_three_maps(arr1, arr2, arr3, titles, suptitle, output_fp, cmap='YlGn'):
    fig, axes = plt.subplots(1, 3, figsize=(18,6), constrained_layout=True)

    ims = []
    for ax, arr, title in zip(axes, [arr1, arr2, arr3], titles):
        im = ax.imshow(arr, cmap=cmap)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.7)
        ims.append(im)

    plt.suptitle(suptitle)
    plt.savefig(output_fp, dpi=300)
    plt.close()

def plot_lai_grid(ds, output_fp, cmap='YlGn'):
    fig, axes = plt.subplots(4, 3, figsize=(16,12), constrained_layout=True)

    months = ["January", "February", "March", "April", "May", "June",
              "July", "August", "September", "October", "November", "December"]

    bands = [f"Band{i+1}" for i in range(12)]

    for i, ax in enumerate(axes.flat):
        if i < len(bands):
            arr = ds[bands[i]].values
            im = ax.imshow(arr, cmap=cmap)
            ax.set_title(months[i])
            ax.axis('off')

    fig.colorbar(im, ax=axes, orientation='horizontal', fraction=0.04, pad=0.08)
    plt.suptitle("Monthly LAI Maps", fontsize=16)
    plt.savefig(output_fp, dpi=300)
    plt.close()

# ==== Main Execution ====

# Load first band (January) for visualization
arr_clipped = load_raster_band(clipped_lai_fp, band=1)
arr_filled = load_raster_band(clipped_filled_lai_fp, band=1)
arr_reproj = load_raster_band(reprojected_lai_fp, band=1)

# Plot clipped, filled, and reprojected LAI
plot_three_maps(
    arr_clipped,
    arr_filled,
    arr_reproj,
    titles=["Clipped LAI", "Clipped & Filled LAI", "Reprojected LAI"],
    suptitle="LAI Preprocessing Steps",
    output_fp=os.path.join(output_folder, "lai_preprocessing_steps.png")
)

# Load NetCDF LAI
ds = xr.open_dataset(lai_nc_fp)

# Plot monthly LAI maps
plot_lai_grid(
    ds,
    output_fp=os.path.join(output_folder, "lai_monthly_grid.png")
)

print("✅ LAI QC visualizations generated.")