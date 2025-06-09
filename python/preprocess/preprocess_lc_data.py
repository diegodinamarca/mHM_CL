import os
import json
import numpy as np
import geopandas as gpd
import rioxarray
from rasterio.enums import Resampling

# Load config
with open("preprocess_config.json") as f:
    config = json.load(f)

# Paths from config
raster_path = config["land_cover_file"]
roi_path = config["roi_file"]
morph_folder = config["morph_folder"]
lc_folder = config["lc_folder"]
dem_file = os.path.join(morph_folder, "dem.asc")
target_resolution = config["cellsize"]
temp_dir = os.path.join(lc_folder, "temp")
os.makedirs(temp_dir, exist_ok=True)

# Reclassification dictionary
reclass_dict = {
    100: 2, 110: 2, 120: 2, 130: 2, 140: 2, 150: 2,
    200: 1, 210: 1, 211: 1, 212: 1, 220: 1, 221: 1, 222: 1,
    230: 1, 231: 1, 232: 1, 240: 1, 241: 1, 242: 1, 250: 1, 251: 1, 252: 1,
    300: 2, 310: 2, 311: 2, 312: 2, 320: 2, 330: 2,
    400: 2, 410: 2, 420: 1, 430: 2, 440: 2, 450: 2,
    500: 2, 510: 2, 520: 2, 530: 2,
    600: 2, 610: 2, 620: 2, 630: 2, 640: 2,
    800: 3,
    900: 2, 910: 3, 920: 2, 930: 3, 931: 3, 932: 3,
    1010: 3, 1020: 3,
    1200: 3
}

# Load ROI
roi = gpd.read_file(roi_path)

# Step 1: Load LC raster and clip to ROI
lc = rioxarray.open_rasterio(raster_path, masked=True).squeeze()
lc = lc.rio.write_crs("EPSG:4326", inplace=True) if lc.rio.crs is None else lc
lc_clipped = lc.rio.clip(roi.geometry.values, roi.crs, drop=True)

# Save clipped version
clipped_path = os.path.join(temp_dir, "landcover_clipped.tif")
lc_clipped.rio.to_raster(clipped_path)

# Step 2: Reclassify using dict
data = lc_clipped.data
reclassified = np.copy(data)
for old, new in reclass_dict.items():
    reclassified[data == old] = new
lc_clipped.data = reclassified

# Save reclassified version
reclass_path = os.path.join(temp_dir, "landcover_reclassified.tif")
lc_clipped.rio.to_raster(reclass_path)

# Step 3: Reproject to EPSG:4326 (safeguard, should already be EPSG:4326)
lc_clipped = lc_clipped.rio.reproject("EPSG:4326", resampling=Resampling.nearest)

# Step 4: Read DEM and align LC to its grid
dem = rioxarray.open_rasterio(dem_file, masked=True).squeeze()
lc_aligned = lc_clipped.rio.reproject_match(dem, resampling=Resampling.nearest)

# Step 5: Export final LC raster as ASCII Grid
final_lc_path = os.path.join(lc_folder, "landcover.asc")
lc_aligned.rio.to_raster(final_lc_path, driver="AAIGrid", nodata=-9999)

print(f"✅ Final landcover.asc saved to: {final_lc_path}")
print(f"🗂️ Clipped raster saved to: {clipped_path}")
print(f"🗂️ Reclassified raster saved to: {reclass_path}")
