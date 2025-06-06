#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 11:56:02 2025

@author: mhm
"""

import os
import f90nml
import json
import shutil

def update_evaluation_gauges(domain_path):
    config_path = os.path.join(domain_path, "preprocess_config.json")
    with open(config_path, "r") as f:
        config = json.load(f)

    exe_folder = os.path.join(domain_path, config["exe_folder"])
    gauges_folder = os.path.join(domain_path, config["gauge_folder"])
    temp_folder = os.path.join(exe_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)

    nml_path = os.path.join(exe_folder, 'mhm.nml')  # Replace with your .nml filename
    
    backup_path = os.path.join(temp_folder, os.path.basename(nml_path))
    shutil.copy2(nml_path, backup_path)
    print(f"Backup of the original .nml file created at: {backup_path}")
    
    
    # Read the existing namelist
    nml = f90nml.read(nml_path)

    # Retrieve the list of .day files in the gauges directory
    gauge_files = [f for f in os.listdir(gauges_folder) if f.endswith('.day')]
    gauge_files.sort()  # Optional: sort the files for consistent ordering

    # Extract gauge IDs from filenames
    gauge_ids = [int(os.path.splitext(f)[0]) for f in gauge_files]

    # Update the evaluation_gauges section
    nml['evaluation_gauges']['ngaugestotal'] = len(gauge_ids)
    nml['evaluation_gauges']['nogauges_domain'] = [len(gauge_ids)]
    nml['evaluation_gauges']['gauge_id'] = gauge_ids
    nml['evaluation_gauges']['gauge_filename'] = gauge_files

    # Write the updated namelist back to the file
    nml.write(nml_path, force=True)
