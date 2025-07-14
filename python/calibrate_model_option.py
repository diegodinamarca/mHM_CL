#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 11:14:13 2025

@author: mhm
"""
import os
import yaml
import shutil
import f90nml

def update_calibrate_option(domain_path, config_name="preprocess_config.yaml"):
    """
    Rewrites the &mainconfig_mhm_mrm section in mhm.nml using values from preprocess_config.yaml,
    including a dynamic 'optimize' value based on the 'calibrate_model' key in the config file.
    """

    # Paths
    config_path = os.path.join(domain_path, config_name)
    exe_path = os.path.join(domain_path, 'exe')
    nml_path = os.path.join(exe_path, 'mhm.nml')
    temp_dir = os.path.join(exe_path, 'temp')
    backup_path = os.path.join(temp_dir, 'mhm_backup.nml')

    os.makedirs(temp_dir, exist_ok=True)

    # Read config file
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading JSON config: {e}")
        return

    # Read and back up NML
    try:
        nml = f90nml.read(nml_path)
        shutil.copy2(nml_path, backup_path)
        print(f"Backup created at: {backup_path}")
    except Exception as e:
        print(f"Error reading or backing up NML: {e}")
        return

    # Parse and coerce 'calibrate_model'
    optimize_flag = config.get("calibrate_model", False)
    if isinstance(optimize_flag, str):
        optimize_flag = optimize_flag.strip().lower() == "true"
    elif not isinstance(optimize_flag, bool):
        optimize_flag = False

    print("From JSON: calibrate_model =", config.get("calibrate_model"))
    print("Used optimize =", optimize_flag)


    # Define new section (use Python True/False for booleans!)
    new_mainconfig_mhm_mrm = {
        "mhm_file_restartin(1)": "test_domain/restart/",
        "mrm_file_restartin(1)": "test_domain/restart/",
        "resolution_routing(1)": 0.03125,
        "timestep": 1,
        "read_restart": False,
        "optimize": optimize_flag,
        "optimize_restart": False,
        "opti_method": 1,
        "opti_function": 2
    }

    # Replace the section
    nml["mainconfig_mhm_mrm"] = new_mainconfig_mhm_mrm

    # Write updated namelist
    try:
        nml.write(nml_path, force=True)
        print(f"&mainconfig_mhm_mrm section updated in: {nml_path}")
    except Exception as e:
        print(f"Error writing updated NML: {e}")

