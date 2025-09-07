import numpy as np
import os

def write_comparison_results(file1_path, file2_path, diff, output_file):
    """
    Write comparison results to a file.
    
    Args:
        file1_path (str): Path to first .dat file
        file2_path (str): Path to second .dat file
        diff (numpy.ndarray): Full differences array between files
        output_file (str): Path to output file
    """
    print(f"Paths: {file1_path}, {file2_path} \n")
    with open(output_file, 'a') as f:
        f.write(f"File 1: {file1_path}\n")
        f.write(f"File 2: {file2_path}\n")
        f.write("Differences (first 100 rows):\n")
        np.savetxt(f, diff[:100], fmt='%.6e')
        f.write("\n\n")

def compare_dat_files(file1_path, file2_path, tolerance=1e-7, output_to_file=None):
    """
    Compare the first two columns of two .dat files (ignoring sigma),
    ignoring comment lines (starting with #).
    
    Args:
        file1_path (str): Path to first .dat file
        file2_path (str): Path to second .dat file  
        tolerance (float): Tolerance for numerical comparison
    
    Returns:
        bool: True if files match within tolerance, False otherwise
    """
    print(f"Comparing files {file1_path}, {file2_path}")
    # Read data from both files, skipping comment lines
    data1 = np.loadtxt(file1_path, comments='#')
    data2 = np.loadtxt(file2_path, comments='#')
    
    # Check if first two columns match in shape
    if data1.shape[0] != data2.shape[0]:
        print(f"Row count mismatch: {file1_path} has {data1.shape[0]} rows, {file2_path} has {data2.shape[0]} rows")
        return False
    
    # Compare only first two columns with tolerance
    diff = np.abs(data1 - data2)
    mismatch_indices = np.where(diff > tolerance)
    
    if output_to_file:
        write_comparison_results(file1_path, file2_path, diff, output_to_file)
    
    if len(mismatch_indices[0]) > 0:
        # Print first mismatch location
        row, col = mismatch_indices[0][0], mismatch_indices[1][0]
        print(f"Mismatch at data row {row}, column {col}:")
        print(f"  {file1_path}: {data1[row, col]}")
        print(f"  {file2_path}: {data2[row, col]}")
        print(f"  Difference: {diff[row, col]} (tolerance: {tolerance})\n")
        return False
    
    print("Files match within tolerance\n")
    return True

compare_dat_files(
os.path.join("atT/OneD_integrated_WAXS_01/Reduction", 'Hor_all_WAXS.dat'),
r"larger_test/1D/WAXS/Averaged/XHor_scan_Run4X_WAXS_averaged.dat"
)
compare_dat_files(
os.path.join("atT/OneD_integrated_SAXS_01/Reduction", 'Hor_all_SAXS.dat'),
r"larger_test/1D/SAXS/Averaged/XHor_scan_Run4X_SAXS_averaged.dat"
)

compare_dat_files(
os.path.join("atT/OneD_integrated_SAXS_01/Reduction", "Run10_AcOH_RampT22_ctr1_all_SAXS.dat"),
"larger_test/1D/SAXS/Averaged/XRun10_AcOHXT22Xctr1X_SAXS_averaged.dat",
output_to_file="differences.txt"
)

compare_dat_files(
os.path.join("atT/OneD_integrated_WAXS_01/Reduction", "Run10_AcOH_RampT22_ctr1_all_WAXS.dat"),
r"larger_test/1D/WAXS/Averaged/XRun10_AcOHXT22Xctr1X_WAXS_averaged.dat",
output_to_file="differences.txt"
)