# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 10:57:38 2025

@author: rappe
"""

import rasterio
import rasterio.mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os
from rasterio.warp import Resampling
from scipy.ndimage import generic_filter
import matplotlib.pyplot as plt
import tempfile
import shutil
import json


# === CONFIGURACIÓN ===
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)
    
soil_folder = config["soil_folder"]
morph_folder = config["morph_folder"]
header_folder = config["header_folder"]
roi_file = config["roi_file"]
cellsize = config["cellsize"]
remove_temp = config["remove_temp"] # elimina el contenido de temp/ al final

horizons = ["0-5", "5-15", "15-30", "30-60", "60-100", "100-200"]
nodata_fill = -9999
fill_iterations = 20         # número de veces que se aplica el filtro de relleno
window_size = 3             # tamaño de la ventana (nxn)

# === FUNCIONES ===
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

def save_raster_tif(path, array, profile, nodata_value):
    profile_updated = profile.copy()
    profile_updated.update({
        'driver': 'GTiff',
        'dtype': rasterio.float32,
        'count': 1,
        'nodata': nodata_value,
        'compress': 'lzw'
    })
    with rasterio.open(path, 'w', **profile_updated) as dst:
        dst.write(array.astype(np.float32), 1)

def save_comparison_plot(original, filled, hz, prop_name, out_folder):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    vmin = np.nanmin(np.where(original == nodata_fill, np.nan, original))
    vmax = np.nanmax(np.where(original == nodata_fill, np.nan, original))

    for ax, data, title in zip(axes, [original, filled], ['Original', 'Rellenado']):
        im = ax.imshow(np.where(data == nodata_fill, np.nan, data), cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f'{prop_name} - Horizonte {hz}cm')
    os.makedirs(out_folder, exist_ok=True)
    plt.savefig(os.path.join(out_folder, f'{prop_name}_{hz}cm.png'), dpi=150)
    plt.close()
    
def read_clip_and_resample_raster(raster_path, shapes, target_cellsize, resampling_method):
    import rasterio
    import rasterio.mask
    import rasterio.warp
    import numpy as np

    with rasterio.open(raster_path) as src:
        # 📌 1. Recorte con shape
        clipped, clipped_transform = rasterio.mask.mask(
            src, shapes, crop=True, nodata=src.nodata
        )

        clipped_array = clipped[0]
        clipped_meta = src.meta.copy()
        clipped_meta.update({
            "height": clipped_array.shape[0],
            "width": clipped_array.shape[1],
            "transform": clipped_transform
        })

        # 📌 2. Calcular nueva resolución
        left, bottom, right, top = rasterio.transform.array_bounds(
            clipped_array.shape[0], clipped_array.shape[1], clipped_transform
        )

        width_deg = right - left
        height_deg = top - bottom

        new_width = max(1, int(np.ceil(width_deg / target_cellsize)))
        new_height = max(1, int(np.ceil(height_deg / target_cellsize)))

        new_transform = rasterio.transform.from_origin(
            west=left,
            north=top,
            xsize=target_cellsize,
            ysize=target_cellsize
        )

        # 📌 3. Crear array de destino
        resampled = np.empty((new_height, new_width), dtype=src.dtypes[0])

        # 📌 4. Reproyectar (resamplear)
        rasterio.warp.reproject(
            source=clipped_array,
            destination=resampled,
            src_transform=clipped_transform,
            src_crs=src.crs,
            dst_transform=new_transform,
            dst_crs=src.crs,
            resampling=resampling_method
        )

        # 📌 5. Perfil actualizado
        new_meta = clipped_meta.copy()
        new_meta.update({
            "height": new_height,
            "width": new_width,
            "transform": new_transform
        })

        return resampled, src.nodata, new_meta

# Leer geometrías del shapefile
gdf = gpd.read_file(roi_file)
shapes = gdf.geometry.values

# Acumuladores globales
global_combinations = []
global_valid_masks = []
global_profiles = []
global_shapes = []
pixel_id_indices = []
pixel_counts = {}

# === PROCESAR CADA HORIZONTE ===
# Crear carpeta temporal para visualizaciones
temp_folder = os.path.join(morph_folder, "temp")
os.makedirs(temp_folder, exist_ok=True)

for idx, hz in enumerate(horizons):
    print(f"🔄 Procesando horizonte {hz} cm")

    # Leer y recortar/resamplear
    bulkd_path = os.path.join(soil_folder, f"Bulkd.{hz}cm.tif")
    clay_path  = os.path.join(soil_folder, f"Clay.{hz}cm.tif")
    sand_path  = os.path.join(soil_folder, f"Sand.{hz}cm.tif")

    bulkd, nodata_b, profile = read_clip_and_resample_raster(bulkd_path, shapes, cellsize, Resampling.bilinear)
    clay, nodata_c, _        = read_clip_and_resample_raster(clay_path, shapes, cellsize, Resampling.bilinear)
    sand, nodata_s, _        = read_clip_and_resample_raster(sand_path, shapes, cellsize, Resampling.bilinear)
    
    
    # === Aplicar relleno focal
    bulkd_filled = fill_nodata_iterative(bulkd, nodata_b, iterations=fill_iterations)
    clay_filled = fill_nodata_iterative(clay, nodata_c, iterations=fill_iterations)
    sand_filled = fill_nodata_iterative(sand, nodata_s, iterations=fill_iterations)
    
    # Guardar ráster original y rellenado
    save_raster_tif(os.path.join(temp_folder, f"Bulkd_{hz}cm_original.tif"), bulkd, profile, nodata_b)
    save_raster_tif(os.path.join(temp_folder, f"Bulkd_{hz}cm_filled.tif"), bulkd_filled, profile, nodata_b)
    
    save_raster_tif(os.path.join(temp_folder, f"Clay_{hz}cm_original.tif"), clay, profile, nodata_c)
    save_raster_tif(os.path.join(temp_folder, f"Clay_{hz}cm_filled.tif"), clay_filled, profile, nodata_c)
    
    save_raster_tif(os.path.join(temp_folder, f"Sand_{hz}cm_original.tif"), sand, profile, nodata_s)
    save_raster_tif(os.path.join(temp_folder, f"Sand_{hz}cm_filled.tif"), sand_filled, profile, nodata_s)
    
    # === Guardar visualización PNG
    save_comparison_plot(bulkd, bulkd_filled, hz, "Bulkd", temp_folder)
    save_comparison_plot(clay, clay_filled, hz, "Clay", temp_folder)
    save_comparison_plot(sand, sand_filled, hz, "Sand", temp_folder)

    # Redondear
    r1 = np.round(bulkd_filled, 2)
    r2 = np.round(clay_filled, 0)
    r3 = np.round(sand_filled, 0)

    # Máscara de valores válidos
    valid_mask = (
        (~np.isnan(r1)) & (r1 != nodata_fill) &
        (~np.isnan(r2)) & (r2 != nodata_fill) &
        (~np.isnan(r3)) & (r3 != nodata_fill)
    )
    # Apilar y guardar combinaciones
    stacked = np.stack([r1, r2, r3], axis=-1)
    valid_pixels = stacked[valid_mask]

    global_combinations.append(valid_pixels)
    global_valid_masks.append(valid_mask)
    global_profiles.append(profile)
    global_shapes.append(r1.shape)

print("📦 Combinando todos los horizontes...")

# Concatenar todas las combinaciones válidas
all_pixels = np.vstack(global_combinations)

# Obtener combinaciones únicas globales
unique_combinations, inverse_indices = np.unique(all_pixels, axis=0, return_inverse=True)

# Contar nSamples por combinación
for idx in inverse_indices:
    pixel_counts[idx] = pixel_counts.get(idx, 0) + 1

# %%
# === GENERAR RASTERS DE SALIDA ===
start = 0
for i, hz in enumerate(horizons):
    horizon_idx = i + 1  # Para el nombre del archivo (1–6)
    print(f"💾 Exportando raster recodificado para horizonte {hz} cm")

    valid_mask = global_valid_masks[i]
    shape = global_shapes[i]
    profile = global_profiles[i]
    n_valid_pixels = np.sum(valid_mask)

    horizon_indices = inverse_indices[start:start + n_valid_pixels]
    start += n_valid_pixels

    output_raster = np.full(shape, nodata_fill, dtype=np.int32)
    output_raster[valid_mask] = horizon_indices
    output_raster[~valid_mask] = nodata_fill
    output_raster = np.where(np.isnan(output_raster), nodata_fill, output_raster)

    # Actualizar profile
    profile.update(
        dtype=rasterio.int32,
        count=1,
        compress='lzw',
        nodata=nodata_fill,
        driver='AAIGrid'  # Para formato .asc
    )

    output_name = f"soil_class_horizon_0{horizon_idx}.asc"
    output_path = os.path.join(morph_folder, output_name)
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(output_raster, 1)
    
    nrows, ncols = output_raster.shape
    transform = profile['transform']
    
    # Extraer coordenadas de la esquina inferior izquierda
    xllcorner = transform[2]
    yllcorner = transform[5] + nrows * transform[4]  # notar que transform[4] es negativo
    
    # Extraer tamaño del píxel
    cellsize = transform[0]  # o abs(transform[4])
    
    # Crear el texto del encabezado
    header_text = (
        f"ncols         {ncols}\n"
        f"nrows         {nrows}\n"
        f"xllcorner     {xllcorner}\n"
        f"yllcorner     {yllcorner}\n"
        f"cellsize      {cellsize}\n"
        f"NODATA_value  {nodata_fill}\n"
    )

    # Guardar en un archivo
    with open(os.path.join(morph_folder,"header.txt"), "w") as f:
        f.write(header_text)
        
    with open(os.path.join(header_folder,"header_morph.txt"), "w") as f:
        f.write(header_text)
        
print("✅ header.txt guardado")     
print("✅ Todos los ráster recodificados exportados.")


# %%
# === EXPORTAR LOOKUP TABLE FORMATEADA ===

lookup_records = []
for i, (bd, clay, sand) in enumerate(unique_combinations):
    n_samples = pixel_counts.get(i, 0)
    lookup_records.append([i, clay, sand, bd, n_samples])

# Crear DataFrame
df_lookup = pd.DataFrame(
    lookup_records,
    columns=["ID", "CLAY[%]", "SAND[%]", "Bd_mu[gcm-3]", "nSamples"]
)

# Eliminar combinaciones con CLAY, SAND y Bd_mu en cero
df_lookup = df_lookup[
    ~((df_lookup["CLAY[%]"] == 0) & (df_lookup["SAND[%]"] == 0) & (df_lookup["Bd_mu[gcm-3]"] == 0))
]

# Formatear texto final
txt_lines = []
txt_lines.append(f"    nSoil_Types           {len(df_lookup)}")
txt_lines.append("     ID    CLAY[%]    SAND[%]    Bd_mu[gcm-3]    nSamples")
# Limpiar posibles NaNs
df_lookup = df_lookup.fillna(0)

# Asegurar tipos
df_lookup["ID"] = df_lookup["ID"].astype(int)
df_lookup["CLAY[%]"] = df_lookup["CLAY[%]"].astype(int)
df_lookup["SAND[%]"] = df_lookup["SAND[%]"].astype(int)
df_lookup["nSamples"] = df_lookup["nSamples"].astype(int)

# Formatear líneas
for _, row in df_lookup.iterrows():
    line = f"{row['ID']:>9}{row['CLAY[%]']:>11}{row['SAND[%]']:>11}{row['Bd_mu[gcm-3]']:>18.2f}{row['nSamples']:>13}"
    txt_lines.append(line)

lookup_txt_path = os.path.join(morph_folder, "soil_classdefinition_iFlag_soilDB_1.txt")
with open(lookup_txt_path, "w") as f:
    f.write("\n".join(txt_lines))

print("✅ Lookup table exportada como 'soil_classdefinition_iFlag_soilDB_1.txt'")

if remove_temp:
    print("🧹 Eliminando carpeta temporal...")
    shutil.rmtree(temp_folder)
