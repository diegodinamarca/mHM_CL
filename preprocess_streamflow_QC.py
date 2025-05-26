
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 2025

@author: Diego Dinamarca

Script to visualize streamflow data for quality control.
Generates monthly and annual mean streamflow plots and reports missing data statistics.
"""

import os
import json
import pandas as pd
import matplotlib.pyplot as plt

# ==== Load Configuration ====

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config_windows.json")

with open(config_path, "r") as f:
    config = json.load(f)

streamflow_data_file = config["streamflow_data_file"]
out_folder = config["out_folder"]

output_folder = os.path.join(out_folder, "Preprocess Streamflow QC")
os.makedirs(output_folder, exist_ok=True)

# ==== Load and Process Data ====

# Load full streamflow CSV
df = pd.read_csv(streamflow_data_file, parse_dates=[0])
df.set_index(df.columns[0], inplace=True)

# ==== Helper Function for Plotting ====

def plot_timeseries(ts, title, ylabel, output_fp):
    plt.figure(figsize=(10,5))
    plt.plot(ts.index, ts.values, linestyle='-')
    plt.title(title)
    plt.xlabel('Date')
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_fp, dpi=300)
    plt.close()

# ==== Loop over each gauge ID ====

for gauge_id in config["gauge_list"]:
    if gauge_id not in df.columns:
        print(f"⚠️  Gauge ID {gauge_id} not found in streamflow data. Skipping.")
        continue

    series_full = df[gauge_id]
    series_valid = series_full.dropna()

    if series_valid.empty:
        print(f"⚠️  No valid data for gauge {gauge_id}. Skipping.")
        continue

    # Define start and end from valid data
    start_date = series_valid.index.min()
    end_date = series_valid.index.max()

    series_cropped = series_full.loc[start_date:end_date]

    # ==== Resample to Monthly and Annual Streamflow ====

    # Monthly total streamflow (sum of daily values)
    monthly_streamflow = series_cropped.resample('ME').mean()

    # Annual mean streamflow (based on monthly sums)
    annual_mean_streamflow = monthly_streamflow.resample('YE').mean()

    # ==== Plotting ====

    # Plot monthly streamflow
    plot_timeseries(
        monthly_streamflow,
        title=f"Monthly Streamflow - Gauge {gauge_id}",
        ylabel="Monthly Streamflow [units]",
        output_fp=os.path.join(output_folder, f"{gauge_id}_monthly_streamflow.png")
    )

    # Plot annual mean streamflow
    plot_timeseries(
        annual_mean_streamflow,
        title=f"Annual Mean Streamflow - Gauge {gauge_id}",
        ylabel="Annual Mean Streamflow [units]",
        output_fp=os.path.join(output_folder, f"{gauge_id}_annual_mean_streamflow.png")
    )

    # ==== Calculate Missing Data Statistics ====

    total_days = (series_cropped.index.max() - series_cropped.index.min()).days + 1
    total_missing = series_cropped.isna().sum()
    percentage_missing = (total_missing / total_days) * 100

    # Save missing data report
    report_lines = [
        f"Gauge ID: {gauge_id}",
        f"Total days in record: {total_days}",
        f"Total missing values: {total_missing}",
        f"Percentage missing: {percentage_missing:.2f}%"
    ]

    report_path = os.path.join(output_folder, f"{gauge_id}_missing_data_report.txt")
    with open(report_path, 'w') as f:
        f.write("\n".join(report_lines))

    print(f"✅ Streamflow QC plots and missing data report generated for Gauge {gauge_id}.")

print("✅ All gauges processed.")