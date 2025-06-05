# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 12:40:15 2025

@author: rappe
"""

import requests
import geopandas as gpd
from shapely.geometry import box, shape
import os
import time
import json
from tqdm import tqdm

# === CONFIGURACIÓN ===
n = 70  # número de zonas latitudinales
output_geo = "C:/Users/rappe/OneDrive/Documentos/FONDECYT_CAMILA/mhm_snow/DATA/SHP/Geologia/geojson_geologia"
os.makedirs(output_geo, exist_ok=True)

base_url = "https://sdngsig.sernageomin.cl/gissdng/rest/services/Geoportal/GeologiaBase/MapServer"

# === Buscar ID de la capa "Geologia" ===
layers_info = requests.get(f"{base_url}?f=json").json()
geologia_layer = next((l for l in layers_info["layers"] if l["name"].strip().lower() == "geologia"), None)

if not geologia_layer:
    print("❌ No se encontró la capa 'Geologia'")
    exit()

layer_id = geologia_layer["id"]
layer_name = geologia_layer["name"]
print(f"✅ Capa encontrada: '{layer_name}' (ID: {layer_id})")

# === Bounding box aproximado de Chile ===
minx, miny = -75.0, -56.0
maxx, maxy = -66.0, -17.0
step = (maxy - miny) / n

all_features = []
features_invalid = []

# === Descargar por zona con barra de progreso ===
print(f"🚧 Dividiendo Chile en {n} zonas. Descargando...")
for i in tqdm(range(n), desc="Zonas"):
    ymin = miny + i * step
    ymax = ymin + step
    geometry_param = f"{minx},{ymin},{maxx},{ymax}"

    query_url = f"{base_url}/{layer_id}/query"
    params = {
        'where': '1=1',
        'geometry': geometry_param,
        'geometryType': 'esriGeometryEnvelope',
        'spatialRel': 'esriSpatialRelIntersects',
        'outFields': '*',
        'returnGeometry': 'true',
        'f': 'json',
        'returnExceededLimitFeatures': 'true',
        'resultRecordCount': 2000,
        'inSR': 4326
    }

    offset = 0
    while True:
        params['resultOffset'] = offset
        response = requests.get(query_url, params=params)
        data = response.json()

        if "features" not in data or not data["features"]:
            break

        all_features.extend(data["features"])
        offset += 2000
        time.sleep(0.5)

print(f"✅ Total de features descargadas: {len(all_features)}")

# === Convertir a GeoDataFrame ===
valid_features = []
features_invalid = []

for f in all_features:
    geom = f.get("geometry")

    # Convertir "rings" de ArcGIS a formato GeoJSON
    if geom and "rings" in geom:
        try:
            converted_geom = {
                "type": "Polygon",
                "coordinates": geom["rings"]
            }
            geometry = shape(converted_geom)
            valid_features.append({
                "attributes": f["attributes"],
                "geometry": geometry
            })
        except Exception as e:
            print(f"❌ Error al convertir una geometría: {e}")
            features_invalid.append(f)
    else:
        features_invalid.append(f)

# Crear GeoDataFrame
attributes = [f["attributes"] for f in valid_features]
geometries = [f["geometry"] for f in valid_features]
gdf = gpd.GeoDataFrame(attributes, geometry=geometries)
gdf.set_crs("EPSG:3857", inplace=True)  # ← CRS original del servidor ArcGIS
gdf = gdf.to_crs("EPSG:4326")           # ← Conversión a lat/lon para salida estándar

# === Guardar como GeoJSON ===
file_path = os.path.join(output_geo, f"{layer_name.replace(' ', '_')}_zonas_{n}.geojson")
gdf.to_file(file_path, driver="GeoJSON")
print(f"✅ GeoJSON guardado en: {file_path}")

# === Guardar features inválidas (sin geometría o con error) ===
if features_invalid:
    error_path = os.path.join(output_geo, f"{layer_name.replace(' ', '_')}_errores.json")
    with open(error_path, "w", encoding="utf-8") as f:
        json.dump(features_invalid, f, ensure_ascii=False, indent=2)
    print(f"⚠️ Se guardaron {len(features_invalid)} features inválidas en: {error_path}")
else:
    print("✅ Todas las geometrías fueron válidas.")
