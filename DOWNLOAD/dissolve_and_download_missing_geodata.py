# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 19:35:36 2025

@author: rappe
"""

import geopandas as gpd
import pandas as pd
from shapely.errors import TopologicalError
import os
import json

# === CONFIGURACIÓN ===
n = 70  # número de zonas latitudinales
output_geo = "C:/Users/rappe/OneDrive/Documentos/FONDECYT_CAMILA/mhm_snow/DATA/SHP/Geologia/geojson_geologia"
os.makedirs(output_geo, exist_ok=True)

layer_name = "Geologia"
geojson_input = os.path.join(output_geo, f"{layer_name}_zonas_{n}.geojson")

# === Leer GeoJSON descargado previamente ===
print(f"📂 Cargando archivo: {geojson_input}")
gdf = gpd.read_file(geojson_input)
gdf.set_crs("EPSG:4326", inplace=True)  # ← Confirmamos que ya está en EPSG:4326

# === Reparar geometrías con buffer(0) ===
print("🔧 Reparando geometrías con buffer(0)...")
geoms_fixed = []
gdf_reparables = []
gdf_no_reparables = []

for geom, attr in zip(gdf.geometry, gdf.drop(columns="geometry").to_dict("records")):
    try:
        fixed = geom.buffer(0)
        if fixed.is_valid:
            geoms_fixed.append(fixed)
            gdf_reparables.append(attr)
        else:
            gdf_no_reparables.append(attr)
    except TopologicalError:
        gdf_no_reparables.append(attr)

# Crear nuevo GeoDataFrame limpio
gdf_clean = gpd.GeoDataFrame(gdf_reparables, geometry=geoms_fixed)
gdf_clean.set_crs("EPSG:4326", inplace=True)
print(f"✅ Geometrías válidas para disolver: {len(gdf_clean)}")
print(f"❌ Geometrías no reparables: {len(gdf_no_reparables)}")

# === Dissolve ===
if "objectid" in gdf_clean.columns:
    print("🔄 Aplicando dissolve por 'objectid'...")
    gdf_dissolved = gdf_clean.dissolve(by="objectid", as_index=False)
else:
    print("⚠️ No se encontró la columna 'objectid'. No se aplicó dissolve.")
    gdf_dissolved = gdf_clean

# === Eliminar duplicados por 'objectid' después del dissolve ===
if "objectid" in gdf_dissolved.columns:
    print("🔄 Eliminando duplicados por 'objectid' tras dissolve...")
    gdf_dissolved = gdf_dissolved.drop_duplicates(subset="objectid", keep="first")
else:
    print("⚠️ No se encontró la columna 'objectid' tras el dissolve.")

# === Cargar features faltantes + forzar el objectid 7093 ===
missing_path = os.path.join(output_geo, "Geologia_faltantes.geojson")
if os.path.exists(missing_path):
    print("📥 Merging con Geologia_faltantes.geojson + objectid 7093...")
    gdf_missing = gpd.read_file(missing_path)
    gdf_missing.set_crs("EPSG:4326", inplace=True)

    # Forzar la descarga del objectid 7093 (si no está ya)
    if 7093 not in gdf_missing["objectid"].values:
        import requests
        from shapely.geometry import shape
        url = "https://sdngsig.sernageomin.cl/gissdng/rest/services/Geoportal/GeologiaBase/MapServer/6/query"
        params = {
            "objectIds": "7093",
            "outFields": "*",
            "returnGeometry": "true",
            "f": "json",
            "returnExceededLimitFeatures": "true"
        }
        response = requests.get(url, params=params).json()
        features = response.get("features", [])
        for f in features:
            geom = f.get("geometry")
            if geom and "rings" in geom:
                converted_geom = {
                    "type": "Polygon",
                    "coordinates": geom["rings"]
                }
                geometry = shape(converted_geom)
                gdf_7093 = gpd.GeoDataFrame([f["attributes"]], geometry=[geometry], crs="EPSG:3857").to_crs("EPSG:4326")
                gdf_missing = gpd.GeoDataFrame(pd.concat([gdf_missing, gdf_7093], ignore_index=True), crs="EPSG:4326")

    # Merge final
    gdf_final = gpd.GeoDataFrame(pd.concat([gdf_dissolved, gdf_missing], ignore_index=True), crs="EPSG:4326")
    gdf_final = gdf_final.drop_duplicates(subset="objectid", keep="first")
else:
    print("⚠️ No se encontró el archivo Geologia_faltantes.geojson, se usará solo el disuelto.")
    gdf_final = gdf_dissolved

# === Guardar como GeoJSON final ===
file_path = os.path.join(output_geo, f"{layer_name}_completo.geojson")
gdf_final.to_file(file_path, driver="GeoJSON")
print(f"✅ GeoJSON final guardado en: {file_path}")

# === Guardar geometrías no reparables ===
if gdf_no_reparables:
    unreparable_path = os.path.join(output_geo, f"{layer_name}_no_reparables.json")
    with open(unreparable_path, "w", encoding="utf-8") as f:
        json.dump(gdf_no_reparables, f, ensure_ascii=False, indent=2)
    print(f"⚠️ Se guardaron {len(gdf_no_reparables)} geometrías no reparables en: {unreparable_path}")
else:
    print("✅ Todas las geometrías fueron reparadas correctamente.")