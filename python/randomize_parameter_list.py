# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:46:43 2025

@author: rappe
"""

mhm_nml = {
    "interception1": {
        "canopyInterceptionFactor": (0.15, 0.4, 0.15, 1, 1)
    },
    "snow1": {
        "snowTreshholdTemperature": (-2.0, 2.0, 1.0, 1, 1),
        "degreeDayFactor_forest": (0.0001, 4.0, 1.5, 1, 1),
        "degreeDayFactor_impervious": (0.0, 1.0, 0.5, 1, 1),
        "degreeDayFactor_pervious": (0.0, 2.0, 0.5, 1, 1),
        "increaseDegreeDayFactorByPrecip": (0.1, 0.9, 0.5, 1, 1),
        "maxDegreeDayFactor_forest": (0.0, 8.0, 3.0, 1, 1),
        "maxDegreeDayFactor_impervious": (0.0, 8.0, 3.5, 1, 1),
        "maxDegreeDayFactor_pervious": (0.0, 8.0, 4.0, 1, 1)
    },
    "soilmoisture1": {
        "orgMatterContent_forest": (0.0, 20.0, 3.4, 1, 1),
        "orgMatterContent_impervious": (0.0, 1.0, 0.1, 1, 1),
        "orgMatterContent_pervious": (0.0, 4.0, 0.6, 1, 1),
        "PTF_lower66_5_constant": (0.6462, 0.9506, 0.76, 1, 1),
        "PTF_lower66_5_clay": (0.0001, 0.0029, 0.0009, 1, 1),
        "PTF_lower66_5_Db": (-0.3727, -0.1871, -0.264, 1, 1),
        "PTF_higher66_5_constant": (0.5358, 1.1232, 0.89, 1, 1),
        "PTF_higher66_5_clay": (-0.0055, 0.0049, -0.001, 1, 1),
        "PTF_higher66_5_Db": (-0.5513, -0.0913, -0.324, 1, 1),
        "PTF_Ks_constant": (-1.2, -0.285, -0.585, 1, 1),
        "PTF_Ks_sand": (0.006, 0.026, 0.0125, 1, 1),
        "PTF_Ks_clay": (0.003, 0.013, 0.0063, 1, 1),
        "PTF_Ks_curveSlope": (60.96, 60.96, 60.96, 0, 1),
        "rootFractionCoefficient_forest": (0.9, 0.999, 0.97, 1, 1),
        "rootFractionCoefficient_impervious": (0.9, 0.95, 0.93, 1, 1),
        "rootFractionCoefficient_pervious": (0.001, 0.09, 0.02, 1, 1),
        "infiltrationShapeFactor": (1.0, 4.0, 1.75, 1, 1)
    },
    "soilmoisture2": {
        "orgMatterContent_forest": (0.0, 20.0, 3.4, 1, 1),
        "orgMatterContent_impervious": (0.0, 1.0, 0.1, 1, 1),
        "orgMatterContent_pervious": (0.0, 4.0, 0.6, 1, 1),
        "PTF_lower66_5_constant": (0.6462, 0.9506, 0.76, 1, 1),
        "PTF_lower66_5_clay": (0.0001, 0.0029, 0.0009, 1, 1),
        "PTF_lower66_5_Db": (-0.3727, -0.1871, -0.264, 1, 1),
        "PTF_higher66_5_constant": (0.5358, 1.1232, 0.89, 1, 1),
        "PTF_higher66_5_clay": (-0.0055, 0.0049, -0.001, 1, 1),
        "PTF_higher66_5_Db": (-0.5513, -0.0913, -0.324, 1, 1),
        "PTF_Ks_constant": (-1.2, -0.285, -0.585, 1, 1),
        "PTF_Ks_sand": (0.006, 0.026, 0.0125, 1, 1),
        "PTF_Ks_clay": (0.003, 0.013, 0.0063, 1, 1),
        "PTF_Ks_curveSlope": (60.96, 60.96, 60.96, 0, 1),
        "rootFractionCoefficient_forest": (0.9, 0.999, 0.97, 1, 1),
        "rootFractionCoefficient_impervious": (0.9, 0.95, 0.93, 1, 1),
        "rootFractionCoefficient_pervious": (0.001, 0.09, 0.02, 1, 1),
        "infiltrationShapeFactor": (1.0, 4.0, 1.75, 1, 1),
        "jarvis_sm_threshold_c1": (0.0, 1.0, 0.5, 1, 1)
    },
    "soilmoisture3": {
        "orgMatterContent_forest": (0.0, 20.0, 3.4, 1, 1),
        "orgMatterContent_impervious": (0.0, 1.0, 0.1, 1, 1),
        "orgMatterContent_pervious": (0.0, 4.0, 0.6, 1, 1),
        "PTF_lower66_5_constant": (0.6462, 0.9506, 0.76, 1, 1),
        "PTF_lower66_5_clay": (0.0001, 0.0029, 0.0009, 1, 1),
        "PTF_lower66_5_Db": (-0.3727, -0.1871, -0.264, 1, 1),
        "PTF_higher66_5_constant": (0.5358, 1.1232, 0.89, 1, 1),
        "PTF_higher66_5_clay": (-0.0055, 0.0049, -0.001, 1, 1),
        "PTF_higher66_5_Db": (-0.5513, -0.0913, -0.324, 1, 1),
        "PTF_Ks_constant": (-1.2, -0.285, -0.585, 1, 1),
        "PTF_Ks_sand": (0.006, 0.026, 0.0125, 1, 1),
        "PTF_Ks_clay": (0.003, 0.013, 0.0063, 1, 1),
        "PTF_Ks_curveSlope": (60.96, 60.96, 60.96, 0, 1),
        "rootFractionCoefficient_forest": (0.97, 0.985, 0.975, 1, 1),
        "rootFractionCoefficient_impervious": (0.97, 0.985, 0.975, 1, 1),
        "rootFractionCoefficient_pervious": (0.97, 0.985, 0.975, 1, 1),
        "infiltrationShapeFactor": (1.0, 4.0, 1.75, 1, 1),
        "rootFractionCoefficient_sand": (0.001, 0.09, 0.09, 1, 1),
        "rootFractionCoefficient_clay": (0.9, 0.999, 0.98, 1, 1),
        "FCmin_glob": (0.1, 0.2, 0.15, 0, 1),
        "FCdelta_glob": (0.1, 0.4, 0.25, 0, 1),
        "jarvis_sm_threshold_c1": (0.0, 1.0, 0.5, 1, 1)
    },
    "soilmoisture4": {
        "orgMatterContent_forest": (0.0, 20.0, 3.4, 1, 1),
        "orgMatterContent_impervious": (0.0, 1.0, 0.1, 1, 1),
        "orgMatterContent_pervious": (0.0, 4.0, 0.6, 1, 1),
        "PTF_lower66_5_constant": (0.6462, 0.9506, 0.76, 1, 1),
        "PTF_lower66_5_clay": (0.0001, 0.0029, 0.0009, 1, 1),
        "PTF_lower66_5_Db": (-0.3727, -0.1871, -0.264, 1, 1),
        "PTF_higher66_5_constant": (0.5358, 1.1232, 0.89, 1, 1),
        "PTF_higher66_5_clay": (-0.0055, 0.0049, -0.001, 1, 1),
        "PTF_higher66_5_Db": (-0.5513, -0.0913, -0.324, 1, 1),
        "PTF_Ks_constant": (-1.2, -0.285, -0.585, 1, 1),
        "PTF_Ks_sand": (0.006, 0.026, 0.0125, 1, 1),
        "PTF_Ks_clay": (0.003, 0.013, 0.0063, 1, 1),
        "PTF_Ks_curveSlope": (60.96, 60.96, 60.96, 0, 1),
        "rootFractionCoefficient_forest": (0.97, 0.985, 0.975, 1, 1),
        "rootFractionCoefficient_impervious": (0.97, 0.985, 0.975, 1, 1),
        "rootFractionCoefficient_pervious": (0.97, 0.985, 0.975, 1, 1),
        "infiltrationShapeFactor": (1.0, 4.0, 1.75, 1, 1),
        "rootFractionCoefficient_sand": (0.001, 0.09, 0.09, 1, 1),
        "rootFractionCoefficient_clay": (0.9, 0.999, 0.98, 1, 1),
        "FCmin_glob": (0.1, 0.2, 0.15, 0, 1),
        "FCdelta_glob": (0.1, 0.4, 0.25, 0, 1)
    },
    "directRunoff1": {
        "imperviousStorageCapacity": (0.0, 5.0, 0.5, 1, 1)
    },
    "PETminus1": {
        "PET_a_forest": (0.3, 1.3, 0.3, 1, 1),
        "PET_a_impervious": (0.3, 1.3, 0.8, 1, 1),
        "PET_a_pervious": (0.3, 1.3, 1.3, 1, 1),
        "PET_b": (0.0, 1.5, 1.5, 1, 1),
        "PET_c": (-2.0, 0.0, -0.7, 1, 1)
    },
    "PET0": {
        "minCorrectionFactorPET": (0.7, 1.3, 0.9, 1, 1),
        "maxCorrectionFactorPET": (0.0, 0.2, 0.1, 1, 1),
        "aspectTresholdPET": (160.0, 200.0, 180.0, 1, 1)
    },
    "PET1": {
        "minCorrectionFactorPET": (0.7, 1.3, 0.93, 1, 1),
        "maxCorrectionFactorPET": (0.0, 0.2, 0.19, 1, 1),
        "aspectTresholdPET": (160.0, 200.0, 171.0, 1, 1),
        "HargreavesSamaniCoeff": (0.0016, 0.003, 0.0023, 1, 1)
    },
    "PET2": {
        "PriestleyTaylorCoeff": (0.75, 1.75, 1.19, 1, 1),
        "PriestleyTaylorLAIcorr": (-0.5, 0.2, 0.058, 1, 1)
    },
    "PET3": {
        "canopyheigth_forest": (15.0, 40.0, 15.0, 1, 1),
        "canopyheigth_impervious": (0.01, 0.5, 0.02, 1, 1),
        "canopyheigth_pervious": (0.1, 5.0, 0.11, 1, 1),
        "displacementheight_coeff": (0.5, 0.85, 0.64, 1, 1),
        "roughnesslength_momentum_coeff": (0.09, 0.16, 0.095, 1, 1),
        "roughnesslength_heat_coeff": (0.07, 0.13, 0.075, 1, 1),
        "stomatal_resistance": (10.0, 200.0, 56.0, 1, 1)
    },
    "interflow1": {
        "interflowStorageCapacityFactor": (75.0, 200.0, 85.0, 1, 1),
        "interflowRecession_slope": (0.0, 10.0, 7.0, 1, 1),
        "fastInterflowRecession_forest": (1.0, 3.0, 1.5, 1, 1),
        "slowInterflowRecession_Ks": (1.0, 30.0, 15.0, 1, 1),
        "exponentSlowInterflow": (0.05, 0.3, 0.125, 1, 1)
    },
    "percolation1": {
        "rechargeCoefficient": (0.0, 50.0, 35.0, 1, 1),
        "rechargeFactor_karstic": (-5.0, 5.0, -1.0, 1, 1),
        "gain_loss_GWreservoir_karstic": (1.0, 1.0, 1.0, 0, 1)
    },
    "routing1": {
        "muskingumTravelTime_constant": (0.31, 0.35, 0.325, 1, 1),
        "muskingumTravelTime_riverLength": (0.07, 0.08, 0.075, 1, 1),
        "muskingumTravelTime_riverSlope": (1.95, 2.1, 2.0, 1, 1),
        "muskingumTravelTime_impervious": (0.09, 0.11, 0.1, 1, 1),
        "muskingumAttenuation_riverSlope": (0.01, 0.5, 0.3, 1, 1)
    },
    "routing2": {
        "streamflow_celerity": (0.1, 15.0, 1.5, 0, 1)
    },
    "routing3": {
        "slope_factor": (0.1, 100.0, 30.0, 0, 1)
    },
    "neutrons1": {
        "Desilets_N0": (300.0, 2000.0, 1500.0, 0, 1),
        "Desilets_LW0": (0.0, 0.2, 0.1783, 0, 1),
        "Desilets_LW1": (0.0, 0.05, 0.0, 0, 1)
    },
    "neutrons2": {
        "COSMIC_N0": (300.0, 2000.0, 1500.0, 0, 1),
        "COSMIC_N1": (0.01, 10.0, 1.0, 0, 1),
        "COSMIC_N2": (0.01, 10.0, 1.0, 0, 1),
        "COSMIC_alpha0": (0.01, 10.0, 1.0, 0, 1),
        "COSMIC_alpha1": (0.01, 10.0, 1.0, 0, 1),
        "COSMIC_L30": (26.56, 424.78, 106.1942, 0, 1),
        "COSMIC_L31": (-118.3, 200.28, 40.9879, 0, 1),
        "COSMIC_LW0": (0.0, 0.2, 0.1783, 0, 1),
        "COSMIC_LW1": (0.0, 0.05, 0.0, 0, 1)
    },
    "geoparameter": {
        "GeoParam": [
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 0, 1),
            (1.0, 1000.0, 100.0, 0, 1),
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 0, 1),
            (1.0, 1000.0, 100.0, 1, 1),
            (1.0, 1000.0, 100.0, 1, 1)
        ]
    }
}

import random
import json
import os

# Step 1: Load exe_folder from config
def get_exe_folder(config_path="preprocess_config_windows.json"):
    """Return the ``exe_folder`` path from a preprocessing configuration file.

    Parameters
    ----------
    config_path : str, optional
        Path to the JSON configuration file. Defaults to
        ``"preprocess_config_windows.json"``.

    Returns
    -------
    str
        The ``exe_folder`` location defined inside the configuration file.
    """

    with open(config_path, "r") as f:
        config = json.load(f)
    return config["exe_folder"]

# Step 2: Randomize value field in parameter tuples
def randomize_mhm_values(nml_dict):
    """Randomize parameter values within the provided namelist dictionary.

    Parameters
    ----------
    nml_dict : dict
        Dictionary representing the ``mhm_parameter.nml`` structure. Tuples of the
        form ``(lb, ub, value, flag, scaling)`` are expected for each parameter.

    Returns
    -------
    None
        The dictionary is modified in place with new random values sampled from
        the lower and upper bounds.
    """

    for block, params in nml_dict.items():
        for key, val in params.items():
            if isinstance(val, tuple) and len(val) == 5:
                lb, ub, _, flag, scaling = val
                new_val = random.uniform(lb, ub)
                params[key] = (lb, ub, new_val, flag, scaling)
            elif key == "GeoParam" and isinstance(val, list):
                for i in range(len(val)):
                    lb, ub, _, flag, scaling = val[i]
                    new_val = random.uniform(lb, ub)
                    val[i] = (lb, ub, new_val, flag, scaling)

# Step 3: Write to NML file
def write_nml_from_dict(nml_dict, output_path):
    """Write the provided namelist dictionary to an ``.nml`` file.

    Parameters
    ----------
    nml_dict : dict
        Dictionary with the same structure as ``mhm_parameter.nml``.
    output_path : str
        Destination path for the generated NML file.

    Returns
    -------
    None
        The function creates or overwrites ``output_path`` with the formatted
        namelist content.
    """

    def fmt(x, width=10, prec=4):
        """Format numbers or strings for consistent column width."""
        return f"{x:>{width}.{prec}f}" if isinstance(x, float) else f"{x:>{width}}"

    with open(output_path, "w") as f:
        for block, params in nml_dict.items():
            f.write(f"&{block}\n")
            if block == "geoparameter" and "GeoParam" in params:
                for i, entry in enumerate(params["GeoParam"], 1):
                    f.write(f"    GeoParam({i:>2},:) =  {fmt(entry[0])}, {fmt(entry[1])}, {fmt(entry[2])}, {fmt(entry[3], 2)}, {fmt(entry[4], 2)}\n")
            else:
                for param, vals in params.items():
                    lb, ub, val, flag, scaling = vals
                    f.write(f"    {param:<35}=  {fmt(lb)}, {fmt(ub)}, {fmt(val)}, {fmt(flag, 2)}, {fmt(scaling, 2)}\n")
            f.write("/\n\n")

# Step 4: Apply and write
exe_folder = get_exe_folder()
os.makedirs(exe_folder, exist_ok=True)
output_file = os.path.join(exe_folder, "mhm_parameter_randomized.nml")

randomize_mhm_values(mhm_nml)
write_nml_from_dict(mhm_nml, output_file)

print(f"✅ Randomized mhm_parameter.nml written to: {output_file}")