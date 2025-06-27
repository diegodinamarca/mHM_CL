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
    """Update the evaluation gauge list inside ``mhm.nml`` for a domain.

    Parameters
    ----------
    domain_path : str
        Path to the domain directory containing ``preprocess_config.json`` and
        gauge files.

    Returns
    -------
    None
        The function modifies ``mhm.nml`` by injecting gauge IDs and filenames
        and saves a backup of the original file.
    """

    config_path = os.path.join(domain_path, "preprocess_config.json")
    with open(config_path, "r") as f:
        config = json.load(f)

    exe_folder = os.path.join(domain_path, config["exe_folder"])
    gauges_folder = os.path.join(domain_path, config["gauge_folder"])
    temp_folder = os.path.join(exe_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)

    nml_path = os.path.join(exe_folder, 'mhm.nml')
    backup_path = os.path.join(temp_folder, os.path.basename(nml_path))
    shutil.copy2(nml_path, backup_path)
    print(f"Backup of the original .nml file created at: {backup_path}")

    # Read the existing namelist
    nml = f90nml.read(nml_path)

    # Retrieve and sort .day gauge files
    gauge_files = sorted([
        f for f in os.listdir(gauges_folder)
        if f.endswith('.day') and not f.startswith("._")
    ])
    gauge_ids = [int(os.path.splitext(f)[0]) for f in gauge_files]

    # Create Fortran-style keys
    gauge_id_keys = {f"gauge_id(1,{i+1})": gauge_ids[i] for i in range(len(gauge_ids))}
    gauge_filename_keys = {f"gauge_filename(1,{i+1})": gauge_files[i] for i in range(len(gauge_files))}

    # Replace the evaluation_gauges section
    nml["evaluation_gauges"] = {
        "ngaugestotal": len(gauge_ids),
        "nogauges_domain(1)": len(gauge_ids),
        **gauge_id_keys,
        **gauge_filename_keys
    }

    # Write the updated namelist
    nml.write(nml_path, force=True)
    print(f"Updated mhm.nml with {len(gauge_ids)} gauges.")
