# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to compare original DEM vs resampled DEM products.
Generates maps, differences, basic histograms, and reports grid size and resolution differences.
"""

import os
import json
import matplotlib.pyplot as plt
import rasterio
from rasterio.mask import mask
import numpy as np
import geopandas as gpd
import pandas as pd

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

original_dem = config["dem_file"]
roi_file = config["roi_file"]
processed_folder = config["morph_folder"]
out_folder = config["out_folder"]
output_folder = os.path.join(out_folder, "Preprocess DEM QC")
os.makedirs(output_folder, exist_ok=True)

variables = {
    "dem": "DEM Elevation (m)",
    "slope": "Slope (degrees)",
    "aspect": "Aspect (degrees)",
    "fdir": "Flow Direction (D8)",
    "facc": "Flow Accumulation (cells)"
}

# ==== Helper functions ====

def load_raster(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        profile = src.profile
        transform = src.transform
    arr = np.where(arr == profile['nodata'], np.nan, arr)
    return arr, profile, transform

def plot_single(arr, title, output_base):
    fig, ax = plt.subplots(figsize=(8,6))
    im = ax.imshow(arr, cmap='terrain')
    fig.colorbar(im, ax=ax, shrink=0.7)
    ax.set_title(title)
    plt.savefig(f'{output_base}.png', dpi=300)
    plt.close()

def plot_difference(orig, resamp, title, output_base):
    diff = resamp - orig

    fig, axes = plt.subplots(2, 1, figsize=(10,12), constrained_layout=True)

    # Difference Map
    im0 = axes[0].imshow(diff, cmap='RdBu', vmin=-np.nanmax(abs(diff)), vmax=np.nanmax(abs(diff)))
    axes[0].set_title(f'Difference (Resampled - Original) {title}')
    fig.colorbar(im0, ax=axes[0], shrink=0.7)

    # Histogram
    axes[1].hist(diff.flatten(), bins=50, color='gray', edgecolor='black')
    axes[1].set_title(f'Difference Histogram {title}')
    axes[1].set_xlabel('Difference Value')
    axes[1].set_ylabel('Frequency')

    plt.savefig(f'{output_base}_diff.png', dpi=300)
    plt.close()

# ==== Main Execution ====

# Clip original DEM to ROI
roi = gpd.read_file(roi_file)
with rasterio.open(original_dem) as src:
    clipped_array, clipped_transform = mask(src, roi.geometry, crop=True)
    clipped_meta = src.meta.copy()
    clipped_meta.update({
        "height": clipped_array.shape[1],
        "width": clipped_array.shape[2],
        "transform": clipped_transform
    })
    orig_dem_arr = clipped_array[0]
    orig_dem_arr = np.where(orig_dem_arr == src.nodata, np.nan, orig_dem_arr)
    orig_dem_res = (src.res[0], src.res[1])

# Plot original clipped DEM
plot_single(orig_dem_arr, "Original DEM (Clipped to ROI)", os.path.join(output_folder, "original_dem"))

# Grid comparison table
comparison_data = []

# Process each variable
for var, var_title in variables.items():
    proc_fp = os.path.join(processed_folder, f"{var}.asc")
    if not os.path.exists(proc_fp):
        print(f"⚠️  Processed file for {var} not found. Skipping.")
        continue

    resamp_arr, resamp_profile, resamp_transform = load_raster(proc_fp)

    # Plot resampled variable
    plot_single(resamp_arr, f"Resampled {var_title}", os.path.join(output_folder, f"{var}_resamp"))

    # Collect grid info
    comparison_data.append({
        "Variable": var,
        "Original shape": f"{orig_dem_arr.shape}",
        "Processed shape": f"{resamp_arr.shape}",
        "Original resolution (deg)": f"{orig_dem_res}",
        "Processed resolution (deg)": f"{(resamp_transform[0], -resamp_transform[4])}"
    })

    # Only compare DEM for differences if shapes match
    if var == "dem":
        if orig_dem_arr.shape == resamp_arr.shape:
            plot_difference(orig_dem_arr, resamp_arr, var_title, os.path.join(output_folder, var))
        else:
            print(f"⚠️ Skipping difference plot for {var} due to different shapes: {orig_dem_arr.shape} vs {resamp_arr.shape}")

    print(f"✅ Visualization completed for {var}.")

# Save comparison table
comparison_df = pd.DataFrame(comparison_data)
comparison_df.to_csv(os.path.join(output_folder, "grid_comparison_summary.csv"), index=False)

print("✅ All comparisons complete.")