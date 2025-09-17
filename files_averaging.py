"""
Run complete file averaging, using fnmatch patterns in variable averaging_patterns, taking those
files, and averaging them together

TODO: Figure out naming conventions for averaged files
REVISIT: Should I average CSVs or just PDIs?


TODO: Store single row of metadata values.
"""
from pathlib import Path
import numpy as np
import fnmatch

def group_files_by_fnmatch_patterns(file_list: list[Path], pattern_groups: list) -> dict:
    """
    Groups files based on fnmatch pattern groups (similar to Step1 notebook approach)
    Input: file_list: List of file paths
           pattern_groups: List of pattern lists (e.g., [["*Run6*RampT*"], ["*Run7*RampT*"]]).
           Each pattern corresponds to a single group
    Output: Dictionary mapping group names to list of file paths
    """
    groups = {}
    matched_files = set()  # Track which files have been matched to avoid duplicates
    
    for i, pattern_group in enumerate(pattern_groups):
        group_files = []
        
        # For each pattern in the group, find matching files
        for pattern in pattern_group:
            for file_path in file_list:
                if file_path not in matched_files and fnmatch.fnmatch(file_path.name, pattern):
                    group_files.append(file_path)
                    matched_files.add(file_path)
        
        # Only create a group if files were found
        if group_files:
            # Create a meaningful group name from the patterns
            group_name = "_".join(pattern_group).replace("*", "X").replace("[", "").replace("]", "")
            groups[group_name] = group_files
        else:
            print(f"Warning: group {group_name} did not match any patterns:")
    unmatched_files = [f for f in file_list if f not in matched_files]
    if unmatched_files:
        print(f"{len(unmatched_files)} files did not match any patterns:")
            
    return groups


def read_dat_file(file_path: Path) -> tuple:
    """
    Reads a .dat file and extracts header comments and data
    Returns: (header_lines, q_values, intensity_values, sigma_values)
    """
    header_lines = []
    data_lines = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.strip())
            else:
                # Skip empty lines
                if line.strip():
                    data_lines.append(line.strip())
    
    # Parse data (skip the column header line if it exists in data_lines)
    data_arrays = []
    for line in data_lines:
        # Skip column header lines that might not start with #
        if 'q_nm' in line or 'sigma' in line:
            continue
        try:
            values = [float(x) for x in line.split()]
            if len(values) == 3:  # q, I, sigma
                data_arrays.append(values)
        except ValueError:
            continue  # Skip malformed lines
    
    if not data_arrays:
        raise ValueError(f"No valid data found in {file_path}")
    
    data_arrays = np.array(data_arrays)
    q_values = data_arrays[:, 0]
    intensity_values = data_arrays[:, 1]
    sigma_values = data_arrays[:, 2]
    
    return header_lines, q_values, intensity_values, sigma_values


def average_dat_files(file_paths: list) -> tuple:
    """
    Averages multiple .dat files using simple column averaging
    Input: List of .dat file paths to average
    Output: (header_lines, avg_q, avg_intensity, avg_sigma)
    """
    if not file_paths:
        raise ValueError("No files to average")
    
    # Use header from first file
    header_lines, first_q, first_intensity, first_sigma = read_dat_file(file_paths[0])
    
    # Initialize arrays for averaging
    all_q = [first_q]
    all_intensity = [first_intensity]
    all_sigma = [first_sigma]
    
    # Read remaining files
    for file_path in file_paths[1:]:
        _, q_vals, intensity_vals, sigma_vals = read_dat_file(file_path)
        all_q.append(q_vals)
        all_intensity.append(intensity_vals)
        all_sigma.append(sigma_vals)
    
    # Convert to numpy arrays and average
    all_q = np.array(all_q)
    all_intensity = np.array(all_intensity)
    all_sigma = np.array(all_sigma)
    
    # Simple averaging across files (axis=0)
    avg_q = np.mean(all_q, axis=0)
    avg_intensity = np.mean(all_intensity, axis=0) 
    avg_sigma = np.mean(all_sigma, axis=0)
    
    return header_lines, avg_q, avg_intensity, avg_sigma


def save_averaged_data(header_lines: list, q, intensity, sigma, output_path: Path):
    """
    Saves averaged data in .dat format with original header
    """
    with open(output_path, 'w') as f:
        # Write header comments from first file
        for line in header_lines:
            f.write(f"{line}\n")
                
        # Write data
        for i in range(len(q)):
            f.write(f"  {q[i]:e}    {intensity[i]:e}   {sigma[i]:e}\n")


def average_files_in_directory(input_dir: Path, detector_type: str, pattern_groups: list):
    """
    Takes in an individual SAXS/Reduction or WAXS/Reduction directory and performs averaging
    for files and writes averaged files to SAXS/Averaged or WAXS/Averaged
    """
    # Get all .dat files in the directory
    dat_files = list(input_dir.glob("*.dat"))
    
    if not dat_files:
        print(f"No .dat files found in {input_dir}")
        return
    
    file_groups = group_files_by_fnmatch_patterns(dat_files, pattern_groups)

    
    # Create output directory
    output_dir = input_dir.parent / "Averaged"
    output_dir.mkdir(exist_ok=True)
    
    print(f"Processing {len(file_groups)} groups in {input_dir}")
    
    # Process each group
    for group_key, files_in_group in file_groups.items():
        print(f"  Averaging {len(files_in_group)} files for pattern: {group_key}")
        
        # Average the files
        header_lines, avg_q, avg_intensity, avg_sigma = average_dat_files(files_in_group)
        
        # Create output filename
        output_filename = f"{group_key}_{detector_type}_averaged.dat"
        output_path = output_dir / output_filename
        
        # Save averaged data
        save_averaged_data(header_lines, avg_q, avg_intensity, avg_sigma, output_path)
        
        print(f"    Saved: {output_path}")
            
def process_directory(base_dir: Path):
    """
    Main function to process both SAXS and WAXS directories using fnmatch patterns
    Input: Path to 1D/ directory containing SAXS/ and WAXS/ subdirectories
    """
    saxs_reduction_dir = base_dir / "SAXS" / "Reduction"
    waxs_reduction_dir = base_dir / "WAXS" / "Reduction"
    
    print(f"Processing directory: {base_dir}")
    print(f"Using averaging patterns: {len(averaging_patterns)} pattern groups")
    for i, patterns in enumerate(averaging_patterns):
        print(f"  Group {i+1}: {patterns}")
    
    if saxs_reduction_dir.exists():
        print("\nProcessing SAXS files...")
        average_files_in_directory(saxs_reduction_dir, "SAXS", averaging_patterns)
    else:
        print(f"SAXS Reduction directory not found: {saxs_reduction_dir}")
    
    if waxs_reduction_dir.exists():
        print("\nProcessing WAXS files...")
        average_files_in_directory(waxs_reduction_dir, "WAXS", averaging_patterns)
    else:
        print(f"WAXS Reduction directory not found: {waxs_reduction_dir}")
    
    print("\nAveraging complete!")

# Process the entire 1D directory

averaging_patterns = [
['*Hor_scan_Run4*'],
['*Run10_AcOH*T20*ctr0*'],
['*Run10_AcOH*T22*ctr1*'],
['*Run11*PS_AcOH*T157*ctr55*']
]
process_directory(Path("larger_test/1D/")
)
