# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 16:26:56 2025

@author: rappe
"""

import json
import os
import datetime

def load_config(json_file):
    """Load configuration from JSON file"""
    with open(json_file, 'r') as f:
        return json.load(f)

def create_mhm_namelist(config):
    """
    Generate mhm.nml from Python data structures
    using output path from JSON config
    """
    # ======================================================================
    # Define all namelist sections as Python dictionaries
    # ======================================================================
    
    namelists = {
        "project_description": {
            "project_details": "mHM model run",
            "setup_description": "Automatically generated configuration",
            "simulation_type": "historical simulation",
            "Conventions": "",
            "contact": "ddinamarcamuller@gmail.com",
            "mHM_details": "Helmholtz Center for Environmental Research - UFZ, Department Computational Hydrosystems, Stochastic Hydrology Group",
            "history": f"Created {datetime.datetime.now().strftime('%Y-%m-%d')}"
        },
        
        "mainconfig": {
            "iFlag_cordinate_sys": 1,
            "nDomains": 2,
            "resolution_Hydrology": [0.03125,0.03125],
            "L0Domain": [1,1],
            "write_restart": True,
            "read_opt_domain_data": [0,0]
        },
        "mainconfig_mhm_mrm": {
            "mhm_file_RestartIn": [f"{config['out_folder']}/mHM_restart_001.nc"],
            "mrm_file_RestartIn": [f"{config['out_folder']}/mRM_restart_001.nc"],
            "resolution_Routing": [0.03125,0.03125],
            "timestep": 1,
            "read_restart": False,
            "mrm_read_river_network": False,
            "optimize": False,
            "optimize_restart": False,
            "opti_method": 1,
            "opti_function": 1
        },
        
        "mainconfig_mrm": {
            "ALMA_convention": True,
            "varnametotalrunoff": "total_runoff",
            "filenametotalrunoff": "total_runoff",
            "gw_coupling": False
        },
        
        "config_riv_temp": {
            "albedo_water": 0.15,
            "pt_a_water": 1.26,
            "emissivity_water": 0.96,
            "turb_heat_ex_coeff": 20.0,
            "max_iter": 50,
            "delta_iter": 1.0e-02,
            "step_iter": 5.0,
            "riv_widths_file": "Q_bkfl",
            "riv_widths_name": "P_bkfl",
            "dir_riv_widths": [
                f"{config['opt_folder']}",
                f"{config['opt_folder']}"
            ]
        },
                
        "directories_general": {
            "dirConfigOut": f"{config['out_folder']}/",
            "dirCommonFiles": f"{config['morph_folder']}/",
            "dir_Morpho": [f"{config['morph_folder']}/"],
            "dir_LCover": [f"{config['lc_folder']}/"],
            "mhm_file_RestartOut": [f"{config['out_folder']}/mHM_restart_001.nc"],
            "mrm_file_RestartOut": [f"{config['out_folder']}/mRM_restart_001.nc"],
            "dir_Out": [f"{config['out_folder']}/"],
            "file_LatLon": [f"{config['latlon_folder']}/latlon.nc"]
        },
        
        "directories_mHM": {
            "inputFormat_meteo_forcings": "nc",
            "bound_error": True,
            "dir_Precipitation": [f"{config['meteo_folder']}/"],
            "dir_Temperature": [f"{config['meteo_folder']}/"],
            "dir_ReferenceET": [f"{config['meteo_folder']}/"],
            "dir_MinTemperature": [f"{config['meteo_folder']}/"],
            "dir_MaxTemperature": [f"{config['meteo_folder']}/"],
            "dir_NetRadiation": [f"{config['meteo_folder']}/"],
            "dir_absVapPressure": [f"{config['meteo_folder']}/"],
            "dir_windspeed": [f"{config['meteo_folder']}/"],
            "dir_Radiation": [f"{config['meteo_folder']}/"],
            "time_step_model_inputs": [0, 0]
        },
        
        "directories_mRM": {
            "dir_Gauges": [
                f"{config['gauge_folder']}/",
                f"{config['gauge_folder']}/"
            ],
            "dir_Total_Runoff": [
                f"{config['out_folder']}/",
                f"{config['out_folder']}/"
            ],
            "dir_Bankfull_Runoff": [
                f"{config['opt_folder']}/",
                f"{config['opt_folder']}/"
            ]
        },
        
        "optional_data": {
            "dir_soil_moisture": [f"{config['opt_folder']}/"],
            "nSoilHorizons_sm_input": 1,
            "timeStep_sm_input": -2,
            "dir_neutrons": [f"{config['opt_folder']}/"],
            "dir_evapotranspiration": [f"{config['opt_folder']}/"],
            "timeStep_et_input": -2,
            "dir_tws": [f"{config['opt_folder']}/"],
            "timeStep_tws_input": -2
        },
        
        "processSelection": {
            "processCase": [
                1,  # Interception
                1,  # Snow
                1,  # Soil moisture
                1,  # Direct runoff
                0,  # PET
                1,  # Interflow
                1,  # Percolation
                1,  # Routing
                1,  # Baseflow
                0,  # Neutrons
                0   # River temp
            ]
        },
                
        "LCover": {
            "nLCoverScene": 1,
            "LCoverYearStart": [1960],
            "LCoverYearEnd": [2021],
            "LCoverfName": ["lancover_final.asc"]
        },
        
        "time_periods": {
            "warming_Days": [365, 180],
            "eval_Per": [
                {
                    "yStart": 1990, "mStart": 1, "dStart": 1,
                    "yEnd": 1993, "mEnd": 12, "dEnd": 31
                },
                {
                    "yStart": 1993, "mStart": 1, "dStart": 1,
                    "yEnd": 1993, "mEnd": 12, "dEnd": 31
                }
            ]
        },
        
        "soildata": {
            "iFlag_soilDB": 0,  # 0: classical, 1: horizon-specific
            "tillageDepth": 200,  # mm
            "nSoilHorizons_mHM": 2,
            "soil_Depth": [200]  # mm (only need n-1 values for iFlag_soilDB=0)
        },
        "LAI_data_information": {
            "timeStep_LAI_input": 0,
            "inputFormat_gridded_LAI": "nc"
        },
        
        "LCover_MPR": {
            "fracSealed_cityArea": 0.6
        },
        
        "directories_MPR": {
            "dir_gridded_LAI": [
                f"{config['exe_folder']}/input/lai/",
                f"{config['exe_folder']}/input/lai/"
            ]
        },
        "evaluation_gauges": {
            "nGaugesTotal": 1,
            "NoGauges_domain": [1],
            "Gauge_id": [[1]],  # 2D array: [domain][gauge]
            "gauge_filename": [["gauge_001.txt"]]
        },
        "inflow_gauges": {
            "nInflowGaugesTotal": 0,
            "NoInflowGauges_domain": [0],
            "InflowGauge_id": [[-9]],
            "InflowGauge_filename": [[""]],
            "InflowGauge_Headwater": [[False]]
        },
        
        "panEvapo": {
            "evap_coeff": [1.30, 1.20, 0.72, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.50]
        },
        
        "nightDayRatio": {
            "read_meteo_weights": False,
            "fnight_prec": [0.46, 0.50, 0.52, 0.51, 0.48, 0.50, 0.49, 0.48, 0.52, 0.56, 0.50, 0.47],
            "fnight_pet": [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
            "fnight_temp": [-0.76, -1.30, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.60, -0.94, -0.53],
            "fnight_ssrd": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "fnight_strd": [0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45]
        },
        
        "Optimization": {
            "nIterations": 7,
            "seed": 1235876,
            "dds_r": 0.2,  # DDS perturbation rate
            "sa_temp": -9.0,  # Simulated Annealing initial temp
            "sce_ngs": 2,  # SCE number of complexes
            "sce_npg": -9,
            "sce_nps": -9,
            "mcmc_opti": False,
            "mcmc_error_params": [0.01, 0.6]
        },
        
        "baseflow_config": {
            "BFI_calc": True
        }
    }

    # ======================================================================
    # Create output directory and write namelist
    # ======================================================================
    os.makedirs(config['exe_folder'], exist_ok=True)
    output_path = os.path.join(config['exe_folder'], 'mhm.nml')

    with open(output_path, 'w') as f:
        # Write header
        f.write(f"! mHM configuration namelist - mhm.nml\n")
        f.write(f"! Complete configuration with all sections\n")
        f.write(f"! Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Write each section
        for section_name, section_data in namelists.items():
            f.write(f"&{section_name}\n")
            
            for key, value in section_data.items():
                if isinstance(value, bool):
                    f.write(f"  {key} = {'.TRUE.' if value else '.FALSE.'}\n")
                elif isinstance(value, list):
                    if all(isinstance(item, dict) for item in value):
                        for i, struct in enumerate(value, 1):
                            for subkey, subvalue in struct.items():
                                f.write(f"  {key}({i})%{subkey} = {subvalue}\n")
                    else:
                        for i, item in enumerate(value, 1):
                            if isinstance(item, list):
                                for j, subitem in enumerate(item, 1):
                                    f.write(f"  {key}({i},{j}) = {subitem}\n")
                            else:
                                f.write(f"  {key}({i}) = {item}\n")
                elif isinstance(value, str):
                    f.write(f'  {key} = "{value}"\n')
                else:
                    f.write(f"  {key} = {value}\n")
            
            f.write("/\n\n")

    print(f"Successfully created complete mhm.nml in: {output_path}")

if __name__ == "__main__":
    # Load configuration
    config = load_config('preprocess_config_windows.json')
    
    # Generate complete namelist
    create_mhm_namelist(config)