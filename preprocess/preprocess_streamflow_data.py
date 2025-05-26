# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 12:26:33 2025

@author: Diego Dinamarca
"""

import pandas as pd
from datetime import datetime
import os
import json

# Load configuration from JSON
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")
with open(config_path, "r") as f:
    config = json.load(f)

streamflow_data_file = config["streamflow_data_file"]
gauge_list = config["gauge_list"]
gauge_folder = config["gauge_folder"]


# Read CSV assuming the first column is the date
df = pd.read_csv(streamflow_data_file, parse_dates=[0])
df.set_index(df.columns[0], inplace=True)

for i, gauge_id in enumerate(gauge_list, start=1):
    if gauge_id not in df.columns:
        print(f"Gauge ID {gauge_id} not found in data.")
        continue

    series_full = df[gauge_id]
    series_valid = series_full.dropna()

    if series_valid.empty:
        print(f"No valid data for gauge {gauge_id}. Skipping.")
        continue

    start_date = series_valid.index.min()
    end_date = series_valid.index.max()

    # Create complete daily time index
    full_range = pd.date_range(start=start_date, end=end_date, freq='D')
    filled_series = series_full.reindex(full_range).fillna(-9999)

    # Header lines
    output_lines = [
        f"{gauge_id} Gauge {i} (daily discharge)",
        "nodata   -9999",
        "n       1       measurements per day [1, 1440]",
        f"start  {start_date:%Y %m %d %H %M}",
        f"end    {end_date:%Y %m %d %H %M}"
    ]

    # Data lines
    for date, value in filled_series.items():
        line = f"{date:%Y  %m  %d  %H  %M}   {value:8.3f}"
        output_lines.append(line)

    output_file = os.path.join(gauge_folder, f"{gauge_id}.txt")
    with open(output_file, 'w') as f:
        f.write("\n".join(output_lines))

    print(f"Wrote file: {output_file}")