
# -*- coding: utf-8 -*-
"""
Refactored script using config file, with temp folder for intermediate files.
"""

import os
import glob
import json
import shutil
from tqdm import tqdm
import spatial_utils as sp
from cdo import Cdo
import xarray as xr
import numpy as np
import geopandas as gpd

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)

roi_path = config["roi_file"]
DATE1 = config.get("start_date", "1960-01-01")
DATE2 = config.get("end_date", "1965-12-31")
remove_temp = config.get("remove_temp", True)
variables = config.get("variables_clim", {})
cellsize_clim = config.get("cellsize_clim")
cellunit = config.get("cellunit")
clim_folder = config.get("meteo_folder")
header_folder = config.get("header_folder")

# Prepare folders
os.makedirs(clim_folder, exist_ok=True)
temp_folder = os.path.join(clim_folder, "temp")
os.makedirs(temp_folder, exist_ok=True)

# Read ROI
gdf = gpd.read_file(roi_path)
ex_gdf = sp.create_extent(gdf)
ex_gdf.to_file(os.path.join(temp_folder,"extent.geojson"), driver="GeoJSON")

# Get roi buffer
# Create a 5 km buffer around features in a shapefile
buffered_gdf = sp.buffer_geodata(ex_gdf, buffer_km=6)
# Save the result
buffered_gdf.to_file(os.path.join(temp_folder, "buffered.geojson"), driver="GeoJSON")

# Convert GeoDataFrame to bounding box for slicing
gdf_bounds = buffered_gdf.total_bounds
X1, Y1, X2, Y2 = gdf_bounds
#%%

# Initialize CDO
cdo = Cdo()

for VAR, paths in tqdm(variables.items(), desc='Variables', unit='var'):
    input_dir = paths.get("input_dir")
    if not input_dir or not clim_folder:
        continue  # Skip if input/output not defined

    print(f"🔄 Processing variable: {VAR}")



    input_files = glob.glob(os.path.join(input_dir, "*.nc"))
    print(f"  📂 Found {len(input_files)} files in {input_dir}")

    processed_files = []
    roi_files = []
    
    # Convert GeoDataFrame to bounding box for slicing
    gdf_bounds = buffered_gdf.total_bounds
    X1, Y1, X2, Y2 = gdf_bounds
    
    # Crop and select variable
    for input_nc in tqdm(input_files, desc=f'Processing {VAR}', unit='file'):
        base_name = os.path.basename(input_nc)
        roi_nc = os.path.join(temp_folder, f"roi_{base_name}")
        var_nc = os.path.join(temp_folder, f"var_{base_name}")

        cdo.sellonlatbox(X1, X2, Y1, Y2, input=input_nc, output=roi_nc)
        cdo.selname(VAR, input=roi_nc, output=var_nc, options="-b F64")

        roi_files.append(roi_nc)
        processed_files.append(var_nc)

    # Merge in time
    if len(processed_files) > 1:
        merged_nc = os.path.join(temp_folder, f"merged_{VAR}.nc")
        cdo.mergetime(input=" ".join(processed_files), output=merged_nc)
    else:
        merged_nc = processed_files[0]

    # Select date range and standardize reference time
    period_nc = os.path.join(temp_folder, f"period_{VAR}.nc")
    final_nc = os.path.join(temp_folder, f"final_{VAR}.nc")
    try:
        cdo.seldate(DATE1, DATE2, input=merged_nc, output=period_nc)
        cdo.setreftime(DATE1, "00:00:00", "days", input=period_nc, output=final_nc)

        # Rename for mHM compatibility
        var_rename_map = {
            "tmean": "tavg",
            "pr": "pre",
            "et0": "pet"
        }
        mhm_var_name = var_rename_map.get(VAR, VAR)
        final_renamed_nc = os.path.join(temp_folder, f"{mhm_var_name}_original_res.nc")
        cdo.chname(f"{VAR},{mhm_var_name}", input=final_nc, output=final_renamed_nc)
        print(f"📛 Renamed variable from '{VAR}' to '{mhm_var_name}' → {final_renamed_nc}")

    except Exception as e:
        print(f"⚠️ No data for date range {DATE1} to {DATE2}. Keeping merged file: {merged_nc}")
        # Clean up if needed
        for f in roi_files + processed_files + [period_nc]:
            if os.path.exists(f):
                os.remove(f)
                print(f"🗑️ Deleted temp file: {f}")
        final_renamed_nc = merged_nc

    # Spatial resampling
    ds = xr.open_dataset(final_renamed_nc)
    lat_min = float(ds.lat.min().values)
    lat_max = float(ds.lat.max().values)
    lon_min = float(ds.lon.min().values)
    lon_max = float(ds.lon.max().values)
    new_lat = np.arange(lat_min, lat_max + cellsize_clim * 0.1, cellsize_clim)
    new_lon = np.arange(lon_min, lon_max + cellsize_clim * 0.1, cellsize_clim)
    new_lat = new_lat[new_lat <= lat_max]
    new_lon = new_lon[new_lon <= lon_max]
    ds_resampled = ds.interp(lat=new_lat, lon=new_lon, method="linear")
    ds_resampled = ds_resampled.fillna(-9999)

    resampled_filename = f"{mhm_var_name}_resampled.nc"
    resampled_path = os.path.join(temp_folder, resampled_filename)
    ds_resampled.to_netcdf(resampled_path)
    print(f"✅ Resampled to resolution: {cellsize_clim} → {resampled_path}")
    
    # Convert GeoDataFrame to bounding box for slicing
    buffered_3km = sp.buffer_geodata(ex_gdf, buffer_km=3)

    ex_bounds = buffered_3km.total_bounds
    X1, Y1, X2, Y2 = ex_bounds
    
    # Copy resampled file without suffix
    simple_filename = f"{mhm_var_name}.nc"
    simple_path = os.path.join(clim_folder, simple_filename)
    cdo.sellonlatbox(X1, X2, Y1, Y2, input=resampled_path, output=simple_path)

    # shutil.copyfile(resampled_path, simple_path)
    print(f"📄 Copied resampled file to: {simple_path}")

    # Write header.txt
    ncols = ds_resampled.sizes["lon"]
    nrows = ds_resampled.sizes["lat"]
    xllcorner = float(ds_resampled.lon.min().values)
    yllcorner = float(ds_resampled.lat.min().values)
    nodata_value = -9999

    header_text = (
        f"ncols         {ncols}\n"
        f"nrows         {nrows}\n"
        f"xllcorner     {xllcorner}\n"
        f"yllcorner     {yllcorner}\n"
        f"cellsize      {cellsize_clim}\n"
        f"NODATA_value  {nodata_value}\n"
    )
    header_path = os.path.join(clim_folder, "header.txt")
    with open(header_path, "w") as f:
        f.write(header_text)
    
    header_path_2 = os.path.join(header_folder, "header_meteo.txt")
    with open(header_path_2, "w") as f:
        f.write(header_text)
    
    print(f"✅ Header saved: {header_path}")

    # Clean up temp folder if requested
    if remove_temp:
        shutil.rmtree(temp_folder)
        print(f"🗑️ Removed all temporary files in: {temp_folder}")