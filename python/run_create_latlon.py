# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 11:57:46 2025

@author: Diego Dinamarca
"""

# run_create_latlon.py
import os
import json
import subprocess

def run_create_latlon(domain_path):
    config_path = os.path.join(domain_path, "preprocess_config.json")

    with open(config_path, "r") as f:
        config = json.load(f)

    header_dir = os.path.join(domain_path, config["header_folder"])
    script_dir = config["latlon_script_folder"]
    output_nc  = os.path.join(domain_path, config["latlon_folder"], "latlon.nc")

    cmd = [
        "python", "create_latlon.py",
        "-c", "epsg:4326",
        "-f", os.path.join(header_dir, "header_morph.txt"),
        "-g", os.path.join(header_dir, "header_meteo.txt"),
        "-e", os.path.join(header_dir, "header_meteo.txt"),
        "-o", output_nc
    ]

    subprocess.run(cmd, cwd=script_dir)