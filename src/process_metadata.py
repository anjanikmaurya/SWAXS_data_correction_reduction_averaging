import re
from pathlib import Path
import pandas as pd


def find_row_number_to_read(raw_file_path: Path):
    """
    Finds the correct row to read metadata from for any given file.
    """
    doc_number_regex = r"_(\d\d\d\d)\."
    capture = re.search(doc_number_regex, str(raw_file_path))
    if capture is None:
        raise RuntimeError(f"Could not find metadata file for raw file {raw_file_path}")
    return int(capture.group(1)) # Represents corresponding index of the CSV file
    
def process_csv_metadata(raw_file_path: Path|str):
    """Takes in a raw data file path, reads a CSV and outputs a tuple containing i0, bstop, and metadata file path"""
    raw_file_path = Path(raw_file_path) # Convert to Path object if not already path
    raw_path = raw_file_path
    raw_filename = raw_path.name 
    csv_dir = raw_path.parent.parent
    
    csv_file = None
    for csv_path in csv_dir.glob("*.csv"):
        if csv_path.with_suffix("").name in raw_filename:
            csv_file = csv_path
            break
    
    if csv_file is None:
        raise RuntimeError(f"process_csv_metadata: Could not find CSV file for {raw_path} in {csv_dir}")
    
    df = pd.read_csv(csv_file)
    index_number_csv = find_row_number_to_read(raw_file_path)
    temp = df.iloc[index_number_csv].to_dict()
    # Keys often have spaces preceding them. Get rid of these.
    proper_dict = {key.strip(): value for key, value in temp.items()}
    return proper_dict


import re
from pathlib import Path
# FIXME: PDI doesn't work properly
def get_saxs_pdi_from_waxs(pdi_path) -> str:
    """Gets SAXS PDI from WAXS PDI. Also would work in vice-versa, but not used"""

    corresponding_pdi_regex = r"(Run\d{1,3}.*scan\d_\d\d\d\d)"
    match = re.search(corresponding_pdi_regex, pdi_path)
    
    if match is None: 
        raise RuntimeError("")
    shared_expression = match.group(1) 
    # Convert to pathlib to change directories
    search_directory = Path(pdi_path).parent.parent / "SAXS"
    for pdi_file in search_directory.glob("*.pdi"):
        if shared_expression in pdi_file.name:
            return str(pdi_file)
    # No file found
    raise FileNotFoundError(f"""Empty WAXS File: {shared_expression} and no 
                    corresponding SAXS found in directory {search_directory}""")
    

def process_pdi_full(raw_file_path: Path, detector_type: str):
    pdi_file_path = str(raw_file_path) + ".pdi"
    
    if not Path(pdi_file_path).exists():
        raise FileNotFoundError(f"PDI file not found: {pdi_file_path}") 
    with open(pdi_file_path, 'r') as f:
        data = f.read()

    if (detector_type.upper() != "SAXS") and (detector_type.upper() != "WAXS"):
        raise RuntimeError("Detector Type must be either SAXS or WAXS")
    if ('All Counters' not in data) and (detector_type.upper() == "WAXS"):
    # Empty PDI file found in WAXS -- check SAXS to see if corresponding file exists
        pdi_file_path = get_saxs_pdi_from_waxs(pdi_file_path)
    counters, motors, extras = get_meta_from_pdi(pdi_file_path)
    if counters is None:
        raise RuntimeError(f"Error Parsing PDI File: {pdi_file_path} -- counters not found")
    return counters['i0'], counters['bstop'], pdi_file_path
        
def get_meta_from_pdi(pdi_file: str):
    # TODO: Fix this function
    """Get motor and counter names and values from PDI file
    Function returns empty dictionary for counters if it cannot process successfully
    """
    with open(pdi_file, 'r') as f:
        data = f.read()
    data = data.replace('\n', ';')
    try:
        counters = re.search('All Counters;(.*);;# All Motors', data).group(1)
        cts = re.split(';|=', counters)
        Counters = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
        motors = re.search('All Motors;(.*);#', data).group(1)
        cts = re.split(';|=', motors)
        Motors = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
    except AttributeError:
        ss1 = '# Diffractometer Motor Positions for image;# '
        ss2 = ';# Calculated Detector Calibration Parameters for image:'
        motors = re.search(f'{ss1}(.*){ss2}', data).group(1)
        cts = re.split(';|=', motors)
        Motors = {c.split()[0]: float(cs) for c, cs in zip(cts[::2], cts[1::2])}
        Motors['TwoTheta'] = Motors['2Theta']
        Counters = {}
    Extras = {}
    if len(data[data.rindex(';') + 1:]) > 0:
        Extras['epoch'] = data[data.rindex(';') + 1:]
    return Counters, Motors, Extras

