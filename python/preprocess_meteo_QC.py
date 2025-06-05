# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 17:15:53 2025

@author: Diego Dinamarca

Extended QC script comparing original vs resampled:
- Long-term spatial means (full & last 30 yr)
- Maps of full-period means
- For pre & pet: monthly -> annual totals, then annual mean
- For tmin, tmax, tavg: daily -> annual mean
Outputs are written to the folder specified by `out_folder` in the config.
"""

import os
import json
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.switch_backend('Agg')  # Disable GUI backend (important for scripts)

# ===== Helper Functions =====

def mask_open(fp, var):
    ds = xr.open_dataset(fp)
    return ds[var].where(ds[var] != -9999)

def long_term_mean(da, start=None, end=None):
    da_sel = da.sel(time=slice(start, end)) if start or end else da
    return da_sel.mean(dim='time')

def full_spatial_mean(da, start=None, end=None):
    return long_term_mean(da, start, end).mean(dim=['lat', 'lon'])

def annual_from_daily(da):
    return da.groupby('time.year').mean(dim='time')

def annual_from_monthly(da):
    mon = da.resample(time='ME').sum()
    return mon.groupby('time.year').mean(dim='time')

def spatial_mean_ts(da_ann):
    return da_ann.mean(dim=['lat', 'lon'])

def plot_map(da2d, title, out):
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    da2d.plot.pcolormesh(
        ax=ax, transform=ccrs.PlateCarree(), cmap='viridis',
        add_colorbar=True, cbar_kwargs={'shrink':0.7}
    )
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title(title)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)

# ===== Load Configuration =====

script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(script_dir, 'preprocess_config_windows.json')) as f:
    cfg = json.load(f)

meteo_folder = cfg['meteo_folder']
out_folder = cfg['out_folder']
qc_folder = os.path.join(out_folder, 'Preprocess meteo QC')
os.makedirs(qc_folder, exist_ok=True)

vars_all = ['pre', 'pet', 'tmin', 'tmax', 'tavg']
period30 = ('1981-01-01', '2021-12-31')

# ===== Main Loop per Variable =====

for var in vars_all:
    # Input file paths
    orig_fp   = os.path.join(meteo_folder, f'{var}_original_res.nc')
    resamp_fp = os.path.join(meteo_folder, f'{var}_resampled.nc')

    da_o = mask_open(orig_fp, var)
    da_r = mask_open(resamp_fp, var)

    # Full and last 30-year spatial means
    full_o = full_spatial_mean(da_o)
    full_r = full_spatial_mean(da_r)
    last30_o = full_spatial_mean(da_o, *period30)
    last30_r = full_spatial_mean(da_r, *period30)

    # Comparison table
    df = pd.DataFrame({
        f'{var}_orig': [float(full_o.values), float(last30_o.values)],
        f'{var}_resamp': [float(full_r.values), float(last30_r.values)],
        'diff': [float((full_r - full_o).values), float((last30_r - last30_o).values)]
    }, index=['all_years','1981_2021'])

    df.round(3).to_csv(os.path.join(qc_folder, f'{var}_longterm_comparison.csv'))

    # Maps of full-period mean
    plot_map(long_term_mean(da_o), f'{var} mean (all) original',
             os.path.join(qc_folder, f'{var}_map_full_orig.png'))
    plot_map(long_term_mean(da_r), f'{var} mean (all) resampled',
             os.path.join(qc_folder, f'{var}_map_full_resamp.png'))

    # Annual time series
    if var in ['pre','pet']:
        ann_o = annual_from_monthly(da_o)
        ann_r = annual_from_monthly(da_r)
    else:
        ann_o = annual_from_daily(da_o)
        ann_r = annual_from_daily(da_r)

    ts_o = spatial_mean_ts(ann_o)
    ts_r = spatial_mean_ts(ann_r)

    df_ts = pd.DataFrame({
        'orig': ts_o.values,
        'resamp': ts_r.values,
        'diff': (ts_r - ts_o).values
    }, index=ts_o['year'].values)

    df_ts.round(3).to_csv(os.path.join(qc_folder, f'{var}_annual_ts_comparison.csv'))

    # Plot time series
    plt.figure(figsize=(8,5))
    plt.plot(ts_o['year'], ts_o, '-o', label='orig')
    plt.plot(ts_r['year'], ts_r, '-x', label='resamp')
    plt.title(f'{var.upper()} annual mean comparison')
    plt.xlabel('Year'); plt.ylabel(var)
    plt.legend(); plt.grid(True)
    plt.savefig(os.path.join(qc_folder, f'{var}_annual_ts.png'), dpi=300, bbox_inches='tight')
    plt.close()
