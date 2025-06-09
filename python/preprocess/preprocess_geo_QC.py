# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025
@author: Diego Dinamarca

Script to visualize geology preprocessing outputs for quality control.
Generates:
- Raster map of reclassified geology classes
- Bar plot of pixel count per geology class
- Vector map of original geology clipped to ROI
"""

import os
import json
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from matplotlib.colors import ListedColormap, BoundaryNorm

# ==== Load Configuration ====
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

out_folder = config["out_folder"]
morph_folder = config["morph_folder"]
original_geo_fp = os.path.join(config["geo_file"])
roi_fp = os.path.join(config["roi_file"])
output_folder = os.path.join(out_folder, "Preprocess Geology QC")
os.makedirs(output_folder, exist_ok=True)

# ==== File Paths ====
geology_asc_fp = os.path.join(morph_folder, "geology_class.asc")
lookup_table_fp = os.path.join(morph_folder, "temp", "LUT_geo_chile.csv")

# ==== Load Data ====
def load_ascii_grid(fp):
    with rasterio.open(fp) as src:
        arr = src.read(1)
        arr = np.where(arr == src.nodata, np.nan, arr)
        profile = src.profile
    return arr, profile

geology_arr, geology_profile = load_ascii_grid(geology_asc_fp)
lookup_df = pd.read_csv(lookup_table_fp)

# ==== Harmonized Colormap ====
n_classes = lookup_df.shape[0]
colors = plt.cm.tab20(np.linspace(0, 1, n_classes))
cmap = ListedColormap(colors)
bounds = np.arange(0.5, n_classes + 1.5)
norm = BoundaryNorm(bounds, cmap.N)
class_labels = lookup_df.iloc[:, 0].tolist()
class_color_dict = dict(zip(class_labels, colors))

# ==== Visualization Functions ====
def plot_geology_comparison(class_arr, cmap, norm, lookup_df, output_fp):
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(class_arr, cmap=cmap, norm=norm)
    ax.set_title("Reclassified Geology Classes")
    cbar = fig.colorbar(im, ax=ax, ticks=np.arange(1, n_classes + 1), shrink=0.7)
    cbar.ax.set_yticklabels(lookup_df.iloc[:, 0].tolist())
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()

def plot_pixel_count_bar(arr, lookup_df, output_fp):
    valid = arr[~np.isnan(arr)].astype(int)
    counts = pd.Series(valid.flatten()).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(counts.index, counts.values, color='skyblue')
    ax.set_xlabel("Reclassified Geology Class")
    ax.set_ylabel("Pixel Count")
    ax.set_title("Pixel Count per Geology Class")
    ax.set_xticks(counts.index)
    ax.set_xticklabels(lookup_df.iloc[:, 0].tolist(), rotation=90)
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()

def plot_original_geology_shapefile(geology_fp, roi_fp, lookup_df, class_color_dict, output_fp):
    gdf = gpd.read_file(geology_fp)
    roi = gpd.read_file(roi_fp)

    if roi.crs != gdf.crs:
        roi = roi.to_crs(gdf.crs)

    gdf_clipped = gpd.clip(gdf, roi)

    lut_col = lookup_df.columns[0]  # e.g., 'GeoClass'
    class_labels = lookup_df[lut_col].tolist()
    shapefile_class_col = 'cd_geol'

    if shapefile_class_col not in gdf_clipped.columns:
        raise ValueError(f"Column '{shapefile_class_col}' not found in original geology shapefile.")

    unmatched = set(gdf_clipped[shapefile_class_col].unique()) - set(class_labels)
    if unmatched:
        print(f"⚠️ Warning: Some geology classes in shapefile not found in lookup table: {unmatched}")

    fig, ax = plt.subplots(figsize=(10, 8))
    for label in class_labels:
        subset = gdf_clipped[gdf_clipped[shapefile_class_col] == label]
        if not subset.empty:
            subset.plot(ax=ax, color=class_color_dict[label], label=str(label), linewidth=0.2, edgecolor='black')

    roi.boundary.plot(ax=ax, color='red', linewidth=1)
    ax.set_title("Original Geology Map (Clipped to ROI)")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize='small')
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()

# ==== Run Visualizations ====
plot_geology_comparison(
    geology_arr,
    cmap,
    norm,
    lookup_df,
    output_fp=os.path.join(output_folder, "geology_class_map.png")
)

plot_pixel_count_bar(
    geology_arr,
    lookup_df,
    output_fp=os.path.join(output_folder, "geology_pixel_count.png")
)

plot_original_geology_shapefile(
    geology_fp=original_geo_fp,
    roi_fp=roi_fp,
    lookup_df=lookup_df,
    class_color_dict=class_color_dict,
    output_fp=os.path.join(output_folder, "original_geology_clipped.png")
)

print("✅ Geology QC visualizations generated successfully.")

