# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to compare original vs resampled soil property rasters.
Generates side-by-side maps.
"""

import os
import json
import rasterio
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

soil_folder = config["soil_folder"]
morph_folder = config["morph_folder"]
roi_file = config["roi_file"]
out_folder = config["out_folder"]

output_folder = os.path.join(out_folder, "Preprocess Soil QC")
os.makedirs(output_folder, exist_ok=True)

horizons = ["0-5", "5-15", "15-30", "30-60", "60-100", "100-200"]
properties = ["Bulkd", "Clay", "Sand"]
nodata_fill = -9999

# ==== Helper Functions ====

def load_clip_raster(raster_path, roi_shapes):
    with rasterio.open(raster_path) as src:
        clipped, _ = mask(src, roi_shapes, crop=True, nodata=src.nodata)
        array = clipped[0]
        array = np.where(array == src.nodata, np.nan, array)
    return array

def load_raster(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        arr = np.where(arr == src.nodata, np.nan, arr)
    return arr

def plot_side_by_side(orig, resamp, title, output_base):
    fig, axes = plt.subplots(1, 2, figsize=(14,6), constrained_layout=True)

    im0 = axes[0].imshow(orig, cmap='viridis')
    axes[0].set_title('Original')
    fig.colorbar(im0, ax=axes[0], shrink=0.7)

    im1 = axes[1].imshow(resamp, cmap='viridis')
    axes[1].set_title('Resampled')
    fig.colorbar(im1, ax=axes[1], shrink=0.7)

    plt.suptitle(title)
    plt.savefig(f'{output_base}.png', dpi=300)
    plt.close()

# ==== Main Execution ====

# Load ROI
roi = gpd.read_file(roi_file)
roi_shapes = roi.geometry.values

# Process each horizon and property
for hz in horizons:
    for prop in properties:
        print(f"🔍 Comparing {prop} at {hz} cm")

        original_fp = os.path.join(soil_folder, f"{prop}.{hz}cm.tif")
        resamp_fp = os.path.join(morph_folder, "temp", f"{prop}_{hz}cm_filled.tif")

        if not os.path.exists(original_fp) or not os.path.exists(resamp_fp):
            print(f"⚠️  Files missing for {prop} {hz} cm. Skipping.")
            continue

        orig_arr = load_clip_raster(original_fp, roi_shapes)
        resamp_arr = load_raster(resamp_fp)

        # Plot side-by-side comparison
        plot_side_by_side(
            orig_arr,
            resamp_arr,
            title=f"{prop} Comparison {hz} cm",
            output_base=os.path.join(output_folder, f"{prop}_{hz}cm")
        )

print("✅ All soil comparisons complete.")