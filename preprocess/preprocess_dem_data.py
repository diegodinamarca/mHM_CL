# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Mon Mar 24 19:01:53 2025

@author: Diego Dinamarca
"""

import os
import json
import shutil
from glob import glob
import richdem as rd
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.enums import Resampling as ResamplingEnum
from rasterio.features import rasterize
import geopandas as gpd
from shapely.geometry import box
from whitebox.whitebox_tools import WhiteboxTools

# ===============================
# Load configuration
# ===============================
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)

input_dem_file = config["dem_file"]
input_hydronet_file = config["hydro_network_file"]
input_roi_file = config["roi_file"]
morph_folder = config["morph_folder"]
projected_crs = config["projected_crs"]
cellsize = config["cellsize"]
burn_network = config["burn_network"]
remove_temp = config["remove_temp"]
ascii_crs = f"EPSG:{config['latlon_crs']}"

os.makedirs(morph_folder, exist_ok=True)
temp_folder = os.path.join(morph_folder, "temp")
os.makedirs(temp_folder, exist_ok=True)

# ===============================
# Helper Functions
# ===============================

def save_rd_array(array, ref_raster_path, output_path):
    with rasterio.open(ref_raster_path) as src:
        profile = src.profile.copy()
        profile.update(dtype='float32', count=1)
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(array.astype('float32'), 1)

def save_as_ascii(input_path, output_path, target_crs, target_res_deg, resampling_method):
    with rasterio.open(input_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, target_crs, src.width, src.height, *src.bounds, resolution=target_res_deg)

        profile = src.profile.copy()
        profile.update({
            'driver': 'AAIGrid',
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'dtype': 'float32',
            'count': 1,
            'nodata': -9999
        })

        with rasterio.open(output_path, 'w', **profile) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=target_crs,
                resampling=resampling_method
            )

# ===============================
# Process First DEM and AOI
# ===============================

# Load AOI
aoi = gpd.read_file(input_roi_file)
aoi = aoi.to_crs(f"EPSG:{projected_crs}")

# Clip DEM
with rasterio.open(input_dem_file) as dem_src:
    dem_crs = dem_src.crs
    aoi_in_dem_crs = aoi.to_crs(dem_crs)
    clipped_array, clipped_transform = mask(dem_src, aoi_in_dem_crs.geometry, crop=True)
    clipped_meta = dem_src.meta.copy()
    clipped_meta.update({
        "height": clipped_array.shape[1],
        "width": clipped_array.shape[2],
        "transform": clipped_transform
    })

    clipped_dem_path = os.path.join(temp_folder, "clipped_dem.tif")
    with rasterio.open(clipped_dem_path, "w", **clipped_meta) as dst:
        dst.write(clipped_array)

# Reproject to projected CRS
projected_dem_path = os.path.join(temp_folder, "projected_dem.tif")
with rasterio.open(clipped_dem_path) as src:
    transform, width, height = calculate_default_transform(
        src.crs, f"EPSG:{projected_crs}", src.width, src.height, *src.bounds)

    profile = src.profile.copy()
    profile.update({
        'crs': f"EPSG:{projected_crs}",
        'transform': transform,
        'width': width,
        'height': height,
        'nodata': -9999  # ✅ Add this line!
    })

    with rasterio.open(projected_dem_path, "w", **profile) as dst:
        reproject(
            source=rasterio.band(src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=f"EPSG:{projected_crs}",
            resampling=Resampling.bilinear
        )

if (burn_network):
    # Load hydro network and reproject to DEM CRS
    hydronet = gpd.read_file(input_hydronet_file).to_crs(f"EPSG:{projected_crs}")
    
    # Clip hydro network to ROI
    hydronet_clipped = gpd.overlay(hydronet, aoi, how='intersection')
    
    # Optional: buffer to widen channels (e.g., 10 meters)
    hydronet_clipped['geometry'] = hydronet_clipped.buffer(10)
    
    # Rasterize the stream network
    with rasterio.open(projected_dem_path) as src:
        dem_array = src.read(1)
        transform = src.transform
        out_meta = src.meta.copy()
    
        stream_raster = rasterize(
            [(geom, 1) for geom in hydronet_clipped.geometry],
            out_shape=dem_array.shape,
            transform=transform,
            fill=0,
            dtype='uint8'
        )
    
        # Burn depth in meters (e.g., lower by 10m where stream exists)
        burn_depth = 10
        dem_array_burned = dem_array.copy()
        dem_array_burned[stream_raster == 1] -= burn_depth
    
        # Save burned DEM
        burned_dem_path = os.path.join(temp_folder, "projected_dem_burned.tif")
        out_meta.update(dtype='float32', nodata=-9999)
    
        with rasterio.open(burned_dem_path, "w", **out_meta) as dst:
            dst.write(dem_array_burned.astype('float32'), 1)

# ===============================
# Hydrological Analysis
# ===============================
if (burn_network):
    rd_dem = rd.LoadGDAL(burned_dem_path, no_data=-9999)
else:
    rd_dem = rd.LoadGDAL(projected_dem_path, no_data=-9999)

# Fill depressions
filled_dem = rd.FillDepressions(rd_dem, in_place=False)

# Save filled DEM
save_rd_array(filled_dem, projected_dem_path, os.path.join(temp_folder, "dem.tif"))

# Slope
slope = rd.TerrainAttribute(filled_dem, attrib='slope_degrees')
save_rd_array(slope, projected_dem_path, os.path.join(temp_folder, "slope.tif"))

# Aspect
aspect = rd.TerrainAttribute(filled_dem, attrib='aspect')
save_rd_array(aspect, projected_dem_path, os.path.join(temp_folder, "aspect.tif"))


wbt = WhiteboxTools()
wbt.set_working_dir(temp_folder)

wbt.fill_depressions(
    dem='projected_dem.tif',
    output='filled.tif'
)

wbt.d8_pointer(
    dem='filled.tif',
    output='fdir.tif'
)

wbt.d8_flow_accumulation(
    i='filled.tif',
    output='facc.tif',
    out_type='cells'
)
print(wbt.tool_help("d8_flow_accumulation"))
# ===============================
# Resample and Export to ASC
# ===============================
res_deg = cellsize

save_as_ascii(os.path.join(temp_folder, "projected_dem.tif"),        os.path.join(morph_folder, "dem.asc"),   ascii_crs, res_deg, Resampling.bilinear)
save_as_ascii(os.path.join(temp_folder, "slope.tif"),      os.path.join(morph_folder, "slope.asc"), ascii_crs, res_deg, Resampling.bilinear)
save_as_ascii(os.path.join(temp_folder, "aspect.tif"),     os.path.join(morph_folder, "aspect.asc"),ascii_crs, res_deg, Resampling.bilinear)
save_as_ascii(os.path.join(temp_folder, "fdir.tif"),   os.path.join(morph_folder, "fdir.asc"),  ascii_crs, res_deg, Resampling.nearest)
save_as_ascii(os.path.join(temp_folder, "facc.tif"), os.path.join(morph_folder, "facc.asc"),  ascii_crs, res_deg, Resampling.bilinear)

# ===============================
# Clean up temp folder
# ===============================

if remove_temp and os.path.exists(temp_folder):
    shutil.rmtree(temp_folder)
    print("✅ Temp folder cleaned.")

print("✅ All processing complete. ASCII grids saved in output folder.")
