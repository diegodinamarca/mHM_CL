import os
import json
import re

# === Load config to get folders ===
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "preprocess_config.json")

with open(config_path, "r") as f:
    config = json.load(f)

morph_folder = config["morph_folder"]
exe_folder = config["exe_folder"]

geology_def_path = os.path.join(morph_folder, "geology_classdefinition.txt")
output_path = os.path.join(exe_folder, "geoparameter_block.nml")

# === Count number of GeoParam entries ===
with open(geology_def_path, "r") as f:
    lines = f.readlines()

geo_lines = [line for line in lines if re.match(r"\s*\d+\s+\d+", line)]
n_geo = len(geo_lines)

# === Format each GeoParam line ===
lb = "1.000"
ub = "1000.00"
val = "100.0"
flag = "1"
scaling = "1"
line_template = "    GeoParam({i:>2},:) =  {lb:>7},  {ub:>9},  {val:>7},  {flag:>3},  {scaling:>3}"

header = "&geoparameter"
footer = "/"
body = "\n".join(
    line_template.format(i=i+1, lb=lb, ub=ub, val=val, flag=flag, scaling=scaling)
    for i in range(n_geo)
)

# === Write output file ===
os.makedirs(exe_folder, exist_ok=True)
with open(output_path, "w") as f:
    f.write(f"{header}\n{body}\n{footer}\n")

print(f"✅ GeoParam block with {n_geo} entries written to:\n{output_path}")

