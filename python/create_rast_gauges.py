# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 10:44:08 2025

@author: Diego Dinamarca
"""

import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from shapely.geometry import Point
import numpy as np
import yaml
import os
import chardet
import csv

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.yaml")

# Load config
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Parameters from config
csv_path = config["fluv_station_file"]
morph_folder = config["morph_folder"]
target_crs = f"EPSG:{config['latlon_crs']}"
nodata_value = -9999

# Path to fdir.asc
fdir_path = os.path.join(morph_folder, "fdir.asc")

# Read the reference raster (fdir.asc)
with rasterio.open(fdir_path) as ref:
    transform = ref.transform
    width = ref.width
    height = ref.height
    raster_crs = ref.crs

# Read station CSV and convert to GeoDataFrame
# Detect encoding
with open(csv_path, 'rb') as rawfile:
    result = chardet.detect(rawfile.read(10000))
encoding = result['encoding']

# Detect delimiter (using csv.Sniffer)
with open(csv_path, 'r', encoding=encoding) as f:
    sample = f.read(1024)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(sample)
    delimiter = dialect.delimiter

# Load CSV using detected encoding and delimiter
df = pd.read_csv(csv_path, encoding=encoding, sep=delimiter)
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df["LON"], df["LAT"]),
    crs=target_crs
)

# Reproject if necessary
if gdf.crs != raster_crs:
    gdf = gdf.to_crs(raster_crs)

# Rasterize using station ID
shapes = ((geom, id_val) for geom, id_val in zip(gdf.geometry, gdf["ID"]))
raster = rasterize(
    shapes=shapes,
    out_shape=(height, width),
    transform=transform,
    fill=nodata_value,
    dtype=np.int32
)

# Save the raster
output_path = os.path.join(morph_folder, "idgauges.asc")

with rasterio.open(
    output_path,
    "w",
    driver="AAIGrid",
    height=height,
    width=width,
    count=1,
    dtype=raster.dtype,
    crs=raster_crs,
    transform=transform,
    nodata=nodata_value,
) as dst:
    dst.write(raster, 1)

print(f"idgauges.asc written and aligned with fdir.asc at {output_path}")