# -*- coding: utf-8 -*-
"""
Procesamiento de Land Cover:
1. Clip con ROI
2. Reclasificación
3. Reproyección a EPSG:4326
4. Exportación como ASCII Grid (.asc) con nodata = -9999
"""

import os
import rasterio
import rasterio.mask
import geopandas as gpd
import numpy as np
import json
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import array_bounds

# === CONFIGURACIÓN ===
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

raster_path = config["land_cover_file"]
roi_path = config["roi_file"]
lc_folder = config["lc_folder"]
target_resolution = 0.001953125  # en grados (~100m)
temp_dir = os.path.join(lc_folder, 'temp')
os.makedirs(temp_dir, exist_ok=True)

# === DICCIONARIO DE RECLASIFICACIÓN ===
reclass_dict = {
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
vec_reclass = np.vectorize(lambda x: reclass_dict.get(x, 0))  # default 0 para valores desconocidos

# === CARGAR Y CLIP ===
roi = gpd.read_file(roi_path)
with rasterio.open(raster_path) as src:
    roi = roi.to_crs(src.crs)
    clipped, clipped_transform = rasterio.mask.mask(src, roi.geometry, crop=True, filled=True)
    clipped = clipped[0]
    meta = src.meta.copy()
    meta.update({
        'height': clipped.shape[0],
        'width': clipped.shape[1],
        'transform': clipped_transform,
        'dtype': 'uint8',
        'nodata': 0,
        'count': 1
    })

# === RECLASIFICAR ===
reclassified = vec_reclass(clipped)

# === REPROYECTAR A EPSG:4326 ===
dst_crs = "EPSG:4326"
bounds = array_bounds(reclassified.shape[0], reclassified.shape[1], meta['transform'])

dst_transform, dst_width, dst_height = calculate_default_transform(
    src_crs=meta['crs'],
    dst_crs=dst_crs,
    width=reclassified.shape[1],
    height=reclassified.shape[0],
    left=bounds[0], bottom=bounds[1], right=bounds[2], top=bounds[3],
    resolution=target_resolution
)

reprojected = np.zeros((dst_height, dst_width), dtype='uint8')
reproject(
    source=reclassified,
    destination=reprojected,
    src_transform=meta['transform'],
    src_crs=meta['crs'],
    dst_transform=dst_transform,
    dst_crs=dst_crs,
    resampling=Resampling.nearest
)

# === APLICAR MASCARA PARA NODATA ===
valid_mask = reclassified != 0
valid_mask_reprojected = np.zeros_like(reprojected, dtype='uint8')

reproject(
    source=valid_mask.astype('uint8'),
    destination=valid_mask_reprojected,
    src_transform=meta['transform'],
    src_crs=meta['crs'],
    dst_transform=dst_transform,
    dst_crs=dst_crs,
    resampling=Resampling.nearest
)

reprojected_with_nodata = np.full_like(reprojected, fill_value=-9999, dtype='int16')
reprojected_with_nodata[valid_mask_reprojected == 1] = reprojected[valid_mask_reprojected == 1].astype('int16')

# === GUARDAR COMO ASCII GRID ===
output_asc_path = os.path.join(lc_folder, 'landcover.asc')

asc_meta = {
    'driver': 'AAIGrid',
    'height': reprojected.shape[0],
    'width': reprojected.shape[1],
    'count': 1,
    'dtype': 'int16',
    'nodata': -9999,
    'crs': dst_crs,
    'transform': dst_transform
}

with rasterio.open(output_asc_path, 'w', **asc_meta) as dst:
    dst.write(reprojected_with_nodata, 1)

print(f"\n✅ Archivo guardado como .asc en: {output_asc_path}")