
import rasterio
from rasterio.transform import from_bounds
from rasterio.features import rasterize
from rasterio.enums import Resampling
import geopandas as gpd
import numpy as np
import pandas as pd
import os
import json
import shutil

def process_vector_to_raster(geo_file, roi_file, target_res, output_folder, clean_temp=True):
    os.makedirs(output_folder, exist_ok=True)
    temp_dir = os.path.join(output_folder, "temp")
    os.makedirs(temp_dir, exist_ok=True)

    # Leer ROI
    roi = gpd.read_file(roi_file)
    if roi.crs is None:
        roi.set_crs("EPSG:4326", inplace=True)

    # Leer shapefile geológico
    geo = gpd.read_file(geo_file)
    if geo.crs is None:
        geo.set_crs("EPSG:4326", inplace=True)

    # Reproyectar geo al CRS de ROI si es necesario
    geo = geo.to_crs(roi.crs)

    # Recortar geología a la extensión de la ROI
    geo_clipped = gpd.clip(geo, roi)

    # Crear LUT para 'cd_geol' (string → número)
    unique_vals = sorted(geo_clipped["cd_geol"].dropna().unique())
    val_map = {val: i+1 for i, val in enumerate(unique_vals)}
    geo_clipped["raster_val"] = geo_clipped["cd_geol"].map(val_map)

    # Crear grilla para rasterización
    bounds = roi.total_bounds
    minx, miny, maxx, maxy = bounds
    width = int((maxx - minx) / target_res)
    height = int((maxy - miny) / target_res)
    transform = from_bounds(minx, miny, maxx, maxy, width, height)

    # Rasterizar
    shapes = zip(geo_clipped.geometry, geo_clipped["raster_val"])
    raster = rasterize(
        shapes,
        out_shape=(height, width),
        transform=transform,
        fill=-9999,
        dtype="int32"
    )

    # Guardar raster como ASCII Grid
    raster_meta = {
        "driver": "AAIGrid",
        "height": raster.shape[0],
        "width": raster.shape[1],
        "count": 1,
        "dtype": "int32",
        "crs": roi.crs,
        "transform": transform,
        "nodata": -9999
    }
    raster_path = os.path.join(output_folder, "geology_class.asc")
    with rasterio.open(raster_path, "w", **raster_meta) as dst:
        dst.write(raster, 1)

    # Guardar LUT
    lut_df = pd.DataFrame({
        "GeoClass": list(val_map.keys()),
        "Value": list(val_map.values())
    })
    lut_df.to_csv(os.path.join(temp_dir, "LUT_geo_chile.csv"), index=False)

    # Crear archivo geology_classdefinition.txt
    n_units = len(val_map)
    param_lines = [f"nGeo_Formations  {n_units}", ""]
    param_lines.append("GeoParam(i)   ClassUnit     Karstic      Description")
    for i, label in enumerate(unique_vals, 1):
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

    with open(os.path.join(output_folder, "geology_classdefinition.txt"), "w") as f:
        f.write("\n".join(param_lines))

    if clean_temp and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        print("🧹 Carpeta 'temp' eliminada.")

    print("✅ Proceso completado: raster generado desde shapefile.")

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")
# Load config
with open(config_path, "r") as f:
    config = json.load(f)
    
geo_file = config["geo_file"]
roi_file = config["roi_file"]
target_res = config["cellsize"]

morph_folder = config["morph_folder"]

process_vector_to_raster(
    geo_file=geo_file,
    roi_file=roi_file,
    target_res=0.01,  # en grados
    output_folder=morph_folder,
    clean_temp=False
)