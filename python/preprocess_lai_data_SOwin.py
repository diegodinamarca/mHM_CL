# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 17:05:15 2025

@author: Diego Dinamarca
"""
import os
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
import geopandas as gpd
from rasterio.shutil import copy as rio_copy
import numpy as np
from scipy.ndimage import generic_filter

def fill_nodata_iterative(array, nodata_value, iterations=3, window_size=3):
    # Asegurar tipo float y trabajar con NaN
    filled = array.astype(np.float32)
    filled[filled == nodata_value] = np.nan  # Reemplazar nodata por np.nan

    def nanmean_filter(values):
        center = values[len(values) // 2]
        if np.isnan(center):
            vals = values[~np.isnan(values)]
            return np.nanmean(vals) if len(vals) > 0 else np.nan
        else:
            return center

    for _ in range(iterations):
        filled = generic_filter(filled, nanmean_filter, size=window_size, mode='constant', cval=np.nan)
    
    # Reemplazar nan por el valor original de nodata para mantener coherencia
    filled[np.isnan(filled)] = nodata_value
    return filled


# Input paths
lai_path = "C:/Users/rappe/OneDrive/Documentos/FONDECYT_CAMILA/mhm_snow/DATA/RAST/LAI/MODIS_LTM_LAI_Chile_uncompressed.tif"
roi_path = "C:/Users/rappe/OneDrive/Documentos/FONDECYT_CAMILA/mhm_snow/DATA/SHP/Cuencas/cuencas_grupos/grupo6.geojson"
target_crs = f"EPSG:{config['latlon_crs']}"
target_resolution = 0.01  # Adjust as needed

# Main output folder
lai_folder = "C:/Users/rappe/OneDrive/Documentos/FONDECYT_CAMILA/mhm_snow/PROC/lai/"
temp_folder = os.path.join(lai_folder, "temp")
os.makedirs(temp_folder, exist_ok=True)

# Output file paths
clipped_path = os.path.join(temp_folder, "clipped_lai.tif")
reprojected_path = os.path.join(temp_folder, "reprojected_lai.tif")
netcdf_path = os.path.join(lai_folder, "lai.nc")

# Load ROI and reproject to raster CRS
roi = gpd.read_file(roi_path)

with rasterio.open(lai_path) as src:
    roi = roi.to_crs(src.crs)
    geoms = [feature["geometry"] for feature in roi.__geo_interface__["features"]]
    out_image, out_transform = mask(src, geoms, crop=True)
    out_meta = src.meta.copy()
    
    # Define nodata value
    nodata_value = src.nodata if src.nodata is not None else -9999
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform,
        "nodata": nodata_value
    })

# Save the original clipped raster (before filling)
with rasterio.open(clipped_path, "w", **out_meta) as dest:
    dest.write(out_image)

# Fill nodata values
out_image_filled = np.empty_like(out_image)
for i in range(out_image.shape[0]):
    out_image_filled[i] = fill_nodata_iterative(out_image[i], nodata_value)

# Save the filled raster
clipped_filled_path = os.path.join(temp_folder, "clipped_filled_lai.tif")
with rasterio.open(clipped_filled_path, "w", **out_meta) as dest:
    dest.write(out_image_filled)

# Reproject clipped raster to target CRS and resolution
with rasterio.open(clipped_path) as src:
    transform, width, height = calculate_default_transform(
        src.crs, target_crs, src.width, src.height, *src.bounds, resolution=target_resolution)

    kwargs = src.meta.copy()
    kwargs.update({
        "crs": target_crs,
        "transform": transform,
        "width": width,
        "height": height
    })

    with rasterio.open(reprojected_path, "w", **kwargs) as dst:
        for i in range(1, src.count + 1):  # Bands 1 to 12
            reproject(
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=target_crs,
                resampling=Resampling.nearest
            )

# Convert final GeoTIFF to NetCDF
with rasterio.open(reprojected_path) as src:
    profile = src.profile.copy()
    profile.update({
        'driver': 'NetCDF',
        'count': 12,
        'dtype': src.dtypes[0],
        'compress': 'deflate',
        'zlevel': 4,
        'crs': target_crs,
        'tiled': True,
        'blockxsize': 128,
        'blockysize': 128,
        'nodata': src.nodata,
        'transform': src.transform
    })

    # Assign "lai" as the variable name
    profile['variables'] = {'lai': {'z': list(range(1, 13))}}

    with rasterio.open(netcdf_path, 'w', **profile) as dst:
        for i in range(1, 13):
            dst.write(src.read(i), i)

print(f"✅ Done! Final NetCDF written to:\n{netcdf_path}")
