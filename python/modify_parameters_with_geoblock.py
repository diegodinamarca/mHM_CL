# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:36:55 2025

@author: rappe
"""

# update_geoparam_block.py

import os
import re
import json

def update_geoparam_block(domain_path):
    # === Load config ===
    config_path = os.path.join(domain_path, "preprocess_config.json")

    with open(config_path, "r") as f:
        config = json.load(f)

    morph_folder = os.path.join(domain_path, config["morph_folder"])
    exe_folder = os.path.join(domain_path, config["exe_folder"])
    os.makedirs(exe_folder, exist_ok=True)

    temp_folder = os.path.join(exe_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)

    # === File paths ===
    geology_def_path = os.path.join(morph_folder, "geology_classdefinition.txt")
    nml_input_path = os.path.join(exe_folder, "mhm_parameter.nml")
    backup_path = os.path.join(temp_folder, "mhm_parameter_backup.nml")

    # === Count number of geology classes ===
    with open(geology_def_path, "r") as f:
        geo_lines = [line for line in f if re.match(r"\s*\d+\s+\d+", line)]
    n_geo = len(geo_lines)

    # === Create new GeoParam block ===
    geo_block = ["&geoparameter"]
    for i in range(1, n_geo + 1):
        geo_block.append(f"    GeoParam({i:>2},:) =   1.000,   1000.00,   100.0,    1,    1")
    geo_block.append("/")
    new_geo_block = "\n".join(geo_block)

    # === Read original namelist ===
    with open(nml_input_path, "r") as f:
        content = f.read()

    # === Backup original file ===
    with open(backup_path, "w") as f:
        f.write(content)

    # === Replace existing GeoParam block ===
    updated_content = re.sub(r"&geoparameter.*?/", new_geo_block, content, flags=re.DOTALL)

    # === Write updated file ===
    with open(nml_input_path, "w") as f:
        f.write(updated_content)

    print(f"✅ Replaced &geoparameter block with {n_geo} entries.")
    print(f"📄 Original namelist backed up to: {backup_path}")
    print(f"✍️  Updated namelist saved to: {nml_input_path}")
