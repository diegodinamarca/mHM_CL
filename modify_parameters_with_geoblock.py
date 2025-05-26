# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:36:55 2025

@author: rappe
"""

import os
import re
import json

# === Load config ===
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)

morph_folder = config["morph_folder"]
exe_folder = config["exe_folder"]
os.makedirs(exe_folder, exist_ok=True)

# === File paths ===
geology_def_path = os.path.join(morph_folder, "geology_classdefinition.txt")
nml_input_path = os.path.join(exe_folder, "mhm_parameter.nml")
nml_output_path = os.path.join(exe_folder, "mhm_parameter_updated.nml")

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

# === Read original namelist and replace block ===
with open(nml_input_path, "r") as f:
    content = f.read()

updated_content = re.sub(
    r"&geoparameter.*?/", new_geo_block, content, flags=re.DOTALL
)

# === Write updated file ===
with open(nml_output_path, "w") as f:
    f.write(updated_content)

print(f"✅ Replaced &geoparameter block with {n_geo} entries.")
print(f"📄 Saved updated namelist to: {nml_output_path}")