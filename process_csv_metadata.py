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
    
def process_csv_metadata(raw_file_path: Path):
    """Takes in a raw data file path, reads a CSV and outputs a tuple containing i0, bstop, and metadata file path"""
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

    i0 = df.iloc[index_number_csv, 3]
    bstop = df.iloc[index_number_csv, 6]
    return i0, bstop, csv_file
