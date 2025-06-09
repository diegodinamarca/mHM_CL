# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 15:38:45 2025

@author: rappe
"""

import rasterio
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.crs import CRS
import geopandas as gpd
import numpy as np
import pandas as pd
import os
import json
import shutil  # Para eliminar carpetas

def process_raster(geo_file, roi_file, target_resolution, output_geo, clean_temp=False):
    os.makedirs(output_geo, exist_ok=True)
    temp_dir = os.path.join(output_geo, "temp")
    os.makedirs(temp_dir, exist_ok=True)

    roi = gpd.read_file(roi_file)
    if roi.crs is None:
        print("No CRS found in ROI. Assigning EPSG:4326 manually...")
        roi.set_crs("EPSG:4326", inplace=True)
    else:
        print(f"ROI CRS: {roi.crs}")

    with rasterio.open(geo_file) as src:
        if src.crs is None:
            raise ValueError("El raster no tiene CRS asignado. No se puede continuar.")
        print(f"Raster CRS: {src.crs}")
        roi = roi.to_crs(src.crs)

        # Buffer de 10 km (aprox. en metros) usando reproyección a UTM
        utm_crs = roi.estimate_utm_crs()
        roi_utm = roi.to_crs(utm_crs)
        roi_buffered_utm = roi_utm.buffer(10000)
        roi_buffered = gpd.GeoDataFrame(geometry=roi_buffered_utm, crs=utm_crs).to_crs(src.crs)

        out_image, out_transform = mask(src, roi_buffered.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })
        buffered_raster = os.path.join(temp_dir, "buffered.tif")
        with rasterio.open(buffered_raster, "w", **out_meta) as dest:
            dest.write(out_image)

    with rasterio.open(buffered_raster) as src:
        res_x, res_y = src.res
        print(f"Resolución original: {res_x} x {res_y}")
        print(f"Resolución objetivo: {target_resolution}")
        scale_x = res_x / target_resolution
        scale_y = res_y / target_resolution
        new_width = int(src.width * scale_x)
        new_height = int(src.height * scale_y)
        print(f"Nuevo tamaño: {new_width} x {new_height}")
        if new_width <= 0 or new_height <= 0:
            raise ValueError("❌ El tamaño del raster reescalado es inválido. Revisa la resolución objetivo.")

        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=Resampling.nearest
        )
        new_transform = src.transform * src.transform.scale(
            (src.width / data.shape[-1]),
            (src.height / data.shape[-2])
        )
        resampled_meta = src.meta.copy()
        resampled_meta.update({
            "height": new_height,
            "width": new_width,
            "transform": new_transform
        })

    resampled_raster = os.path.join(temp_dir, "resampled.tif")
    with rasterio.open(resampled_raster, "w", **resampled_meta) as dest:
        dest.write(data)

    with rasterio.open(resampled_raster) as src:
        roi_cut = roi.to_crs(src.crs)
        out_image, out_transform = mask(src, roi_cut.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })
        cropped_raster = os.path.join(temp_dir, "cropped.tif")
        with rasterio.open(cropped_raster, "w", **out_meta) as dest:
            dest.write(out_image)

    with rasterio.open(cropped_raster) as dest:
        band = dest.read(1).astype("float32")
        nodata_val = dest.nodata if dest.nodata is not None else -9999
        band[band == nodata_val] = np.nan
        unique_vals = np.unique(band[~np.isnan(band)])
        val_map = {val: i+1 for i, val in enumerate(unique_vals)}
        reclassified_band = np.copy(band)
        for original, new in val_map.items():
            reclassified_band[band == original] = new
        reclassified_band[np.isnan(reclassified_band)] = -9999
        reclassified_band = reclassified_band.astype(np.int32)

        reclassified_meta = dest.meta.copy()
        reclassified_meta.update({"nodata": -9999})
        reclassified_raster = os.path.join(output_geo, "geology_class.asc")
        with rasterio.open(reclassified_raster, "w", driver="AAIGrid", **reclassified_meta) as final:
            final.write(reclassified_band[np.newaxis, :, :])

    lut_df = pd.DataFrame({
        "OriginalValue": list(val_map.keys()),
        "ReclassifiedValue": list(val_map.values())
    })
    lut_df.to_csv(os.path.join(temp_dir, "lookup_table.csv"), index=False)

    n_units = len(val_map)
    param_lines = [f"nGeo_Formations  {n_units}", ""]
    param_lines.append("GeoParam(i)   ClassUnit     Karstic      Description")
    for i in range(1, n_units + 1):
        param_lines.append(f"{i:>10} {i:>13} {0:>11}      GeoUnit-{i}")
    param_lines.append("!<-END\n")
    param_lines.append("!***********************************")
    param_lines.append("! NOTES")
    param_lines.append("!***********************************")
    param_lines.append("1 = Karstic")
    param_lines.append("0 = Non-karstic\n")
    param_lines.append("IMPORTANT ::")
    param_lines.append("   Ordering has to be according to the ordering in mhm_parameter.nml")
    param_lines.append("   (namelist: geoparameter)")

    with open(os.path.join(output_geo, "geology_classdefinition.txt"), "w") as f:
        f.write("\n".join(param_lines))

    # ✅ Borrar carpeta temp si clean_temp es True
    if clean_temp and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        print("🧹 Carpeta 'temp' eliminada.")

    print("✅ Proceso completado.")

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")
# Load config
with open(config_path, "r") as f:
    config = json.load(f)
    
geo_file = config["geo_file"]
roi_file = config["roi_file"]
target_res = config["cellsize"]
morph_folder = config["morph_folder"]
# Ejemplo de uso:
process_raster(geo_file, roi_file, target_res, morph_folder)

