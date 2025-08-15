import re
from pathlib import Path
import pandas as pd

def process_csv_metadata(raw_file_path: str):
    """Takes in a raw data file path, reads a CSV and outputs a tuple containing i0, bstop, and metadata file path"""
    # Captures the part of the string Hor_scan_..._scan1. 1 Can be replaced
    regex_identify_metadata = r".*?(Hor_scan.*scan\d)"
    # This string must contain the following letters
    string_in_metadata: re.Match[str] | None = re.search(regex_identify_metadata, raw_file_path)
    if string_in_metadata is None:
        raise RuntimeError(f"Could not find metadata file for raw file {raw_file_path}")
    
    # Loop through previous directory, finding csv file containing string_in_metadata
    raw_path = Path(raw_file_path)
    csv_dir = raw_path.parent.parent.parent
    captured_group = string_in_metadata.group(1)
    
    csv_file = None
    for csv_path in csv_dir.glob("*.csv"):
        if captured_group in csv_path.name:
            csv_file = csv_path
            break
    
    if csv_file is None:
        raise RuntimeError(f"Could not find CSV file containing '{captured_group}' in {csv_dir}")
    
    df = pd.read_csv(csv_file)
    if df.shape[0] > 1:
        raise NotImplementedError("Multiple metadata rows in CSV. Not yet implemented")
    # Todo: don't just take first row
    i0 = pd.to_numeric(df.iloc[0, 3], errors='coerce')
    bstop = pd.to_numeric(df.iloc[0, 6], errors='coerce')

    return i0, bstop, csv_file