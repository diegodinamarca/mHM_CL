import os
import yaml
import shutil
import f90nml
from dateutil import parser as date_parser

def update_time_periods(domain_path):
    """
    Updates the &time_periods section in the mhm.nml file based on the start_date and end_date
    specified in the preprocess_config.yaml file located in the domain_path.

    Parameters:
    - domain_path (str): Path to the domain directory containing preprocess_config.yaml and the exe folder.
    """

    # Paths
    config_path = os.path.join(domain_path, 'preprocess_config.yaml')
    exe_path = os.path.join(domain_path, 'exe')
    nml_path = os.path.join(exe_path, 'mhm.nml')
    temp_dir = os.path.join(exe_path, 'temp')
    backup_path = os.path.join(temp_dir, 'mhm_backup.nml')

    os.makedirs(temp_dir, exist_ok=True)

    # Read config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading JSON: {e}")
        return

    try:
        start_date = date_parser.parse(config["start_date"])
        end_date = date_parser.parse(config["end_date"])
    except Exception as e:
        print(f"Error parsing dates: {e}")
        return

    # Read and back up NML
    try:
        nml = f90nml.read(nml_path)
        shutil.copy2(nml_path, backup_path)
        print(f"Backup created at: {backup_path}")
    except Exception as e:
        print(f"Error reading or backing up NML: {e}")
        return

    # Create clean &time_periods section
    new_time_periods = {
        "warming_days(1)": 1100,
        "eval_per(1)%ystart": start_date.year+5,
        "eval_per(1)%mstart": start_date.month,
        "eval_per(1)%dstart": start_date.day,
        "eval_per(1)%yend": end_date.year,
        "eval_per(1)%mend": end_date.month,
        "eval_per(1)%dend": end_date.day,
    }

    # Replace the section completely
    nml["time_periods"] = new_time_periods

    # Write
    try:
        nml.write(nml_path, force=True)
        print(f"Updated mhm.nml written to: {nml_path}")
    except Exception as e:
        print(f"Error writing updated NML: {e}")
