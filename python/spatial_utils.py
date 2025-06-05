# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 10:24:11 2025

@author: Diego Dinamarca
"""

import geopandas as gpd
from shapely.geometry import box
from typing import Union
import rasterio
from rasterio.mask import mask
import xarray as xr
from shapely.geometry import mapping

def get_extent_points(file_path):
    # Load the file (shapefile, GeoJSON, or KML)
    gdf = gpd.read_file(file_path)

    # Get extent (bounding box) [minx, miny, maxx, maxy]
    minx, miny, maxx, maxy = gdf.total_bounds

    # Bounding box corners
    extent = {
        "bottom_left": (minx, miny),
        "bottom_right": (maxx, miny),
        "top_left": (minx, maxy),
        "top_right": (maxx, maxy),
    }

    return extent

def get_extent(file_path):
    # Load the file (shapefile, GeoJSON, or KML)
    gdf = gpd.read_file(file_path)

    # Get extent (bounding box) [minx, miny, maxx, maxy]
    minx, miny, maxx, maxy = gdf.total_bounds

    # Bounding box corners
    extent = {
        "xmin": minx,
        "xmax": maxx,
        "ymin": miny,
        "ymax": maxy
    }

    return extent


def buffer_geodata(
    data: Union[str, gpd.GeoDataFrame],
    buffer_km: float,
    utm_epsg: int = 32718
) -> gpd.GeoDataFrame:
    """
    Buffers a GeoDataFrame or file (shapefile or GeoJSON) by N kilometers.
    Projects to UTM before buffering and reprojects to EPSG:4326 if original was in lat/lon.

    Parameters:
    - data (str or GeoDataFrame): Path to a shapefile/GeoJSON or a GeoDataFrame.
    - buffer_km (float): Buffer distance in kilometers.
    - utm_epsg (int): EPSG code of the UTM zone to use for buffering (default 32718).

    Returns:
    - GeoDataFrame: Buffered geometry in EPSG:4326 (if original data was lat/lon).
    """
    # Load if a file path is given
    if isinstance(data, str):
        gdf = gpd.read_file(data)
    elif isinstance(data, gpd.GeoDataFrame):
        gdf = data.copy()
    else:
        raise ValueError("Input must be a file path or a GeoDataFrame.")

    # Ensure CRS is defined
    if gdf.crs is None:
        raise ValueError("Input GeoDataFrame must have a defined CRS.")

    # Check if the CRS is geographic (lat/lon)
    is_latlon = gdf.crs.to_epsg() in [4326, 4269, 4258]

    # Project to UTM for buffering
    if is_latlon:
        gdf = gdf.to_crs(epsg=utm_epsg)

    # Apply the buffer (buffer_km converted to meters)
    gdf["geometry"] = gdf.buffer(buffer_km * 1000)

    # Reproject back to EPSG:4326 if original was lat/lon
    if is_latlon:
        gdf = gdf.to_crs(epsg=4326)

    return gdf



def create_extent(
    data: Union[str, gpd.GeoDataFrame]
) -> gpd.GeoDataFrame:
    """
    Reads a shapefile/GeoJSON or a GeoDataFrame and returns a rectangular geometry
    matching the extent of the input features.

    Parameters:
    - data (str or GeoDataFrame): Path to input shapefile/GeoJSON or a GeoDataFrame.

    Returns:
    - GeoDataFrame: A GeoDataFrame with a single rectangular geometry covering the extent.
    """
    # Load if file path is given
    if isinstance(data, str):
        gdf = gpd.read_file(data)
    elif isinstance(data, gpd.GeoDataFrame):
        gdf = data.copy()
    else:
        raise ValueError("Input must be a file path or a GeoDataFrame.")

    if gdf.crs is None:
        raise ValueError("Input must have a defined CRS.")

    # Get the extent (bounding box)
    minx, miny, maxx, maxy = gdf.total_bounds
    rectangle_geom = box(minx, miny, maxx, maxy)

    # Return as GeoDataFrame
    rectangle_gdf = gpd.GeoDataFrame(geometry=[rectangle_geom], crs=gdf.crs)

    return rectangle_gdf




def clip_raster(
    gdf: gpd.GeoDataFrame,
    raster_path: str
) -> Union[rasterio.io.DatasetReader, xr.Dataset]:
    """
    Clips a raster file (GeoTIFF, ASCII, NetCDF) using the geometry of a GeoDataFrame.

    Parameters:
    - gdf (GeoDataFrame): GeoDataFrame with clipping geometry.
    - raster_path (str): Path to the raster file (.tif, .asc, .nc).

    Returns:
    - Clipped raster as:
        - rasterio.io.MemoryFile.DatasetReader for raster (.tif, .asc),
        - xarray.Dataset for NetCDF (.nc).
    """
    # Ensure GeoDataFrame has CRS
    if gdf.crs is None:
        raise ValueError("GeoDataFrame must have a CRS.")

    if raster_path.endswith((".tif", ".asc")):
        with rasterio.open(raster_path) as src:
            # Reproject GeoDataFrame to raster CRS if needed
            if gdf.crs != src.crs:
                gdf = gdf.to_crs(src.crs)

            # Get geometries in GeoJSON-like dict format
            geometries = [mapping(geom) for geom in gdf.geometry]

            # Mask the raster with geometries
            out_image, out_transform = mask(src, geometries, crop=True)
            out_meta = src.meta.copy()

            # Update metadata with new shape and transform
            out_meta.update({
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })

            # Write to in-memory dataset
            memfile = rasterio.io.MemoryFile()
            with memfile.open(**out_meta) as dest:
                dest.write(out_image)
            return memfile.open()

    elif raster_path.endswith(".nc"):
        ds = xr.open_dataset(raster_path)

        # Try guessing spatial dims
        lat_dim = next((d for d in ds.dims if "lat" in d.lower()), None)
        lon_dim = next((d for d in ds.dims if "lon" in d.lower()), None)

        if not lat_dim or not lon_dim:
            raise ValueError("Cannot identify lat/lon dimensions in NetCDF file.")

        # Convert GeoDataFrame to bounding box for slicing
        gdf_bounds = gdf.to_crs("EPSG:4326").total_bounds
        minx, miny, maxx, maxy = gdf_bounds

        # Clip using xarray indexing
        clipped_ds = ds.sel(
            {lon_dim: slice(minx, maxx), lat_dim: slice(maxy, miny)}
        )

        return clipped_ds

    else:
        raise ValueError("Unsupported raster format. Supported: .tif, .asc, .nc")
