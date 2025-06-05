# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to visualize geology preprocessing outputs for quality control.
Generates a map comparing original geology and final reclassified geology classes.
Also generates a bar plot of pixel counts per geology class.
"""

import os
import json
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

out_folder = config["out_folder"]
morph_folder = config["morph_folder"]

output_folder = os.path.join(out_folder, "Preprocess Geology QC")
os.makedirs(output_folder, exist_ok=True)

# ==== File Paths ====

geology_asc_fp = os.path.join(morph_folder, "geology_class.asc")
lookup_table_fp = os.path.join(morph_folder, "temp", "LUT_geo_chile.csv")
original_geo_fp = os.path.join(morph_folder, "temp", "cropped.tif")

# ==== Load Data ====

def load_ascii_grid(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        arr = np.where(arr == src.nodata, np.nan, arr)
        profile = src.profile
    return arr, profile

def load_raster(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        arr = np.where(arr == src.nodata, np.nan, arr)
        profile = src.profile
    return arr, profile

# Load geology rasters
geology_arr, geology_profile = load_ascii_grid(geology_asc_fp)
original_arr, original_profile = load_raster(original_geo_fp)

# Load lookup table
lookup_df = pd.read_csv(lookup_table_fp)

# Create colormap for geology classes
n_classes = lookup_df.shape[0]
colors = plt.cm.tab20(np.linspace(0, 1, n_classes))
cmap = ListedColormap(colors)
bounds = np.arange(0.5, n_classes+1.5)
norm = BoundaryNorm(bounds, cmap.N)

# ==== Visualization ====

def plot_geology_comparison(orig_arr, class_arr, cmap, norm, lookup_df, output_fp):
    original_classes = lookup_df.columns[0]

    fig, axes = plt.subplots(1, 2, figsize=(18,8), constrained_layout=True)

    im0 = axes[0].imshow(orig_arr, cmap='terrain')
    axes[0].set_title("Original Geology Map")
    fig.colorbar(im0, ax=axes[0], shrink=0.7)

    im1 = axes[1].imshow(class_arr, cmap=cmap, norm=norm)
    axes[1].set_title("Reclassified Geology Classes")
    cbar = fig.colorbar(im1, ax=axes[1], ticks=np.arange(1, len(lookup_df)+1), shrink=0.7)
    cbar.ax.set_yticklabels(lookup_df[original_classes].tolist())

    plt.suptitle("Geology Map Comparison", fontsize=16)
    plt.savefig(output_fp, dpi=300)
    plt.close()

def plot_pixel_count_bar(arr, lookup_df, output_fp):
    original_classes = lookup_df.columns[0]

    valid = arr[~np.isnan(arr)].astype(int)
    counts = pd.Series(valid.flatten()).value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(10,6))
    ax.bar(counts.index, counts.values, color='skyblue')
    ax.set_xlabel("Reclassified Geology Class")
    ax.set_ylabel("Pixel Count")
    ax.set_title("Pixel Count per Geology Class")
    ax.set_xticks(counts.index)
    ax.set_xticklabels(lookup_df[original_classes].tolist(), rotation=90)
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()


# Plot geology comparison map
plot_geology_comparison(
    original_arr,
    geology_arr,
    cmap,
    norm,
    lookup_df,
    output_fp=os.path.join(output_folder, "geology_comparison_map.png")
)

# Plot pixel count bar chart
plot_pixel_count_bar(
    geology_arr,
    lookup_df,
    output_fp=os.path.join(output_folder, "geology_pixel_count.png")
)

print("✅ Geology QC comparison and pixel count visualizations generated.")
