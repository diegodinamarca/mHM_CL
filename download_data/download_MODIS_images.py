# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:37:02 2025

@author: Diego Dinamarca
"""
# %% import libraries and initialize GEE
import ee
import geemap
import os
import json
import time

# Authenticate and initialize Earth Engine (only needed once)
try:
    ee.Initialize()
except Exception as e:
    ee.Authenticate()
    ee.Initialize()
    
output_folder = "C:/Users/rappe/OneDrive/Documentos/FONDECYT CAMILA/mhm_snow"  # Windows
# output_folder = "/home/yourusername/Documents"  # Linux/macOS

# write a file to custom directory
# outfile = os.path.join(output_folder, "gee_map.html")

# %% explore dates
mod15_lastimg = ee.ImageCollection("MODIS/061/MOD15A2H").sort('system:time_start', False).first()

# Extract the date
latest_date = ee.Date(mod15_lastimg.get('system:time_start')).format('YYYY-MM-dd')

# Print the latest date
print("Latest available date:", latest_date.getInfo())

# %% Load input area and products
# Load the GeoJSON file
aoi_path = "C:/Users/rappe/OneDrive/Documentos/FONDECYT CAMILA/mhm_snow/DATA/SHP/Chile_Polygon.geojson"  # Replace with your file path

# Read the GeoJSON file
with open(aoi_path) as f:
    aoi_data = json.load(f)

# Convert GeoJSON to an Earth Engine Geometry
aoi = ee.Geometry(aoi_data["features"][0]["geometry"])

# %% Calculate Long-term mean monthly LAI values

# Select a date range
start_date = '2000-01-01'
end_date = '2024-12-31'

# Load MOD15A2H (Terra) and MYD15A2H (Aqua) collections
modis_terra = ee.ImageCollection("MODIS/061/MOD15A2H").filterDate(start_date, end_date).select("Lai_500m")
modis_aqua = ee.ImageCollection("MODIS/061/MYD15A2H").filterDate(start_date, end_date).select("Lai_500m")

# Merge both collections
modis_combined = modis_terra.merge(modis_aqua)

# Function to calculate long-term mean for each month
def monthly_mean(month):
    month = ee.Number(month)
    
    # Filter images for the given month
    monthly_images = modis_combined.filter(ee.Filter.calendarRange(month, month, 'month'))
    
    # Compute the mean LAI for that month
    return monthly_images.mean().set("month", month)

# Generate a list of months (1 to 12)
months = ee.List.sequence(1, 12)

# Apply the function to all months and create a single multi-band image
monthly_images = months.map(lambda m: monthly_mean(m))

# Convert the list of images to an ImageCollection
monthly_images = ee.ImageCollection(monthly_images)

# Stack images into a single multi-band image
final_image = monthly_images.toBands()

# Generate band names dynamically based on available months
num_bands = final_image.bandNames().size()
expected_bands = ee.List(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

# Ensure the rename function gets the right number of bands
final_image = final_image.rename(expected_bands.slice(0, num_bands))


# %% Export image

# Export to Google Drive
export_task = ee.batch.Export.image.toDrive(
    image=final_image,
    description="LongTerm_Monthly_LAI",
    folder="LTM_MODIS_LAI",
    fileNamePrefix="MODIS_LongTerm_Monthly_LAI",
    scale=500,
    region=aoi,
    fileFormat="GeoTIFF",
    maxPixels=1e13
)

# Start the export task
export_task.start()
print("Export started. Check Google Drive for the results.")
# %% check task status
# Function to check task status and print errors if the task fails
def check_task_status(task):
    while task.status()['state'] in ['READY', 'RUNNING']:
        print(f"Task Status: {task.status()['state']}...")
        time.sleep(30)  # Wait 30 seconds before checking again

    # Get final task status
    final_status = task.status()
    print(f"Final Task Status: {final_status['state']}")

    # If the task failed, print the error message
    if final_status['state'] == 'FAILED':
        print(f"Error Message: {final_status.get('error_message', 'No error message available.')}")
    elif final_status['state'] == 'COMPLETED':
        print("Task completed successfully!")

check_task_status(export_task)
