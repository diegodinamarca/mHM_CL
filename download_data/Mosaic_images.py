# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 17:59:46 2025

@author: Diego Dinamarca
"""

import rasterio
from rasterio.merge import merge
import glob
import os

def merge_multiband_images(image_paths, output_path):
    # Abrimos cada imagen y las agregamos a una lista
    src_files = [rasterio.open(fp) for fp in image_paths]
    
    # Usamos rasterio.merge.merge para unir las imágenes
    mosaic, out_trans = merge(src_files)
    
    # Actualizamos los metadatos a partir del primer archivo
    out_meta = src_files[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "count": mosaic.shape[0],
        "compress": "lzw"  # or "deflate"
    })
    
    # Escribimos el mosaico en un archivo nuevo
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)
    
    # Cerramos los archivos abiertos
    for src in src_files:
        src.close()

if __name__ == "__main__":
    # Especifica la ruta donde se encuentran tus imágenes
    carpeta_imagenes = "C:/Users/rappe/OneDrive/Documentos/FONDECYT CAMILA/mhm_snow/DATA/RAST/Topography/DEM/DEM_SRTM_CHILE"
    # Buscamos todos los archivos TIFF en la carpeta
    lista_imagenes = glob.glob(os.path.join(carpeta_imagenes, "*.tif"))
    
    # Definimos el nombre y ruta del archivo de salida
    salida = "SRTM_DEM_Chile.tif"
    
    # Llamamos a la función para unir las imágenes
    merge_multiband_images(lista_imagenes, salida)
    
    print("Mosaico creado y guardado en:", salida)
