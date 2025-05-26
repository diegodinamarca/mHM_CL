# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to visualize land cover preprocessing for quality control.
Generates maps of original cropped, resampled, and reclassified land cover rasters.
Also outputs the lookup table for reclassification.
"""

import os
import json
import geopandas as gpd
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

lc_folder = config["lc_folder"]
lc_file = config["land_cover_file"]

out_folder = config["out_folder"]
roi_file = config["roi_file"]

output_folder = os.path.join(out_folder, "Preprocess Land Cover QC")
os.makedirs(output_folder, exist_ok=True)

original_lc_fp = os.path.join(lc_file)
final_lc_fp = os.path.join(lc_folder, "landcover_final.asc")

# ==== Load ROI and Reproject ====
roi = gpd.read_file(roi_file)

with rasterio.open(original_lc_fp) as src:
    raster_crs = src.crs

roi = roi.to_crs(raster_crs)
roi_shapes = roi.geometry.values


# ==== Reclassification Dictionary ====
reclassify_dict = {
    100: 2, 110: 2, 120: 2, 130: 2, 140: 2, 150: 2,
    200: 1, 210: 1, 211: 1, 212: 1, 220: 1, 221: 1, 222: 1,
    230: 1, 231: 1, 232: 1, 240: 1, 241: 1, 242: 1, 250: 1, 251: 1, 252: 1,
    300: 2, 310: 2, 311: 2, 312: 2, 320: 2, 330: 2,
    400: 2, 410: 2, 420: 1, 430: 2, 440: 2, 450: 2,
    500: 2, 510: 2, 520: 2, 530: 2,
    600: 2, 610: 2, 620: 2, 630: 2, 640: 2,
    800: 3,
    900: 2, 910: 2, 920: 2, 930: 2, 931: 3, 932: 3,
    1010: 3, 1020: 3,
    1200: 3
}
category_labels = {1: "Tree", 2: "No Tree", 3: "Impervious"}
category_colors = ['forestgreen', 'khaki', 'gray']
category_cmap = ListedColormap(category_colors)
category_bounds = [0.5, 1.5, 2.5, 3.5]
category_norm = BoundaryNorm(category_bounds, category_cmap.N)

# ==== Helper Functions ====

def load_clip_raster(raster_fp, shapes):
    with rasterio.open(raster_fp) as src:
        clipped, _ = mask(src, shapes, crop=True)
        arr = clipped[0]
        arr = np.where(arr == src.nodata, np.nan, arr)
        profile = src.profile
    return arr, profile

def load_raster(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        arr = np.where(arr == src.nodata, np.nan, arr)
        profile = src.profile
    return arr, profile

def plot_two_maps(arr1, arr2, title1, title2, suptitle, output_fp, cmap='tab20'):
    fig, axes = plt.subplots(1, 2, figsize=(14,6), constrained_layout=True)

    im0 = axes[0].imshow(arr1, cmap=cmap)
    axes[0].set_title(title1)
    fig.colorbar(im0, ax=axes[0], shrink=0.7)

    im1 = axes[1].imshow(arr2, cmap=cmap)
    axes[1].set_title(title2)
    fig.colorbar(im1, ax=axes[1], shrink=0.7)

    plt.suptitle(suptitle)
    plt.savefig(output_fp, dpi=300)
    plt.close()

def plot_categorical_map(arr, title, output_fp):
    fig, ax = plt.subplots(figsize=(8,6))
    im = ax.imshow(arr, cmap=category_cmap, norm=category_norm)
    cbar = fig.colorbar(im, ax=ax, ticks=[1,2,3], shrink=0.7)
    cbar.ax.set_yticklabels([category_labels[v] for v in [1,2,3]])
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()

# ==== Main Execution ====

# Load original clipped
orig_arr, orig_profile = load_clip_raster(original_lc_fp, roi_shapes)

# Load final reclassified
final_arr, final_profile = load_raster(final_lc_fp)

# ==== Visualization ====

# Plot cropped original vs resampled raster
plot_two_maps(
    orig_arr,
    final_arr,
    title1="Original Land Cover (Cropped)",
    title2="Processed Land Cover (Resampled + Reclassified)",
    suptitle="Land Cover Preprocessing Comparison",
    output_fp=os.path.join(output_folder, "lc_comparison.png")
)

# Plot final reclassified categorical map
plot_categorical_map(
    final_arr,
    title="Final Reclassified Land Cover",
    output_fp=os.path.join(output_folder, "lc_final_reclassified.png")
)

# ==== Output Lookup Table ====

lut = pd.DataFrame(list(reclassify_dict.items()), columns=["Original_Class", "Reclassified_Class"])
lut.sort_values("Original_Class", inplace=True)
lut.to_csv(os.path.join(output_folder, "landcover_reclassification_lookup.csv"), index=False)

print("✅ Land cover QC visualizations and lookup table generated.")