# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 11:57:46 2025

@author: Diego Dinamarca
"""

import subprocess
import os
import json

script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)

header_dir = config["header_folder"]
script_dir = config["latlon_script_folder"]
output_nc  = config["latlon_folder"]

cmd = [
    "python", "create_latlon.py",
    "-c", "epsg:4326",
    "-f", os.path.join(header_dir, "header_morph.txt"),
    "-g", os.path.join(header_dir, "header_meteo.txt"),
    "-e", os.path.join(header_dir, "header_meteo.txt"),
    "-o", output_nc
]

subprocess.run(cmd, cwd=script_dir)