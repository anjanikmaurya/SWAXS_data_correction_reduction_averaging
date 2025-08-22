import numpy as np

def compare_dat_files(file1_path, file2_path, tolerance=1e-7):
    """
    Compare two .dat files, ignoring comment lines (starting with #).
    
    Args:
        file1_path (str): Path to first .dat file
        file2_path (str): Path to second .dat file  
        tolerance (float): Tolerance for numerical comparison
    
    Returns:
        bool: True if files match within tolerance, False otherwise
    """
    
    # Read data from both files, skipping comment lines
    data1 = np.loadtxt(file1_path, comments='#')
    data2 = np.loadtxt(file2_path, comments='#')
    
    # Check if shapes match
    if data1.shape != data2.shape:
        print(f"Shape mismatch: {file1_path} has shape {data1.shape}, {file2_path} has shape {data2.shape}")
        return False
    
    # Compare data with tolerance
    diff = np.abs(data1 - data2)
    mismatch_indices = np.where(diff > tolerance)
    
    if len(mismatch_indices[0]) > 0:
        # Print first mismatch location
        row, col = mismatch_indices[0][0], mismatch_indices[1][0]
        print(f"Mismatch at data row {row}, column {col}:")
        print(f"  {file1_path}: {data1[row, col]}")
        print(f"  {file2_path}: {data2[row, col]}")
        print(f"  Difference: {diff[row, col]} (tolerance: {tolerance})")
        return False
    
    print("Files match within tolerance")
    return True

# compare_dat_files('atT_smaller/OneD_integrated_SAXS_01/Reduction/Hor_all_SAXS.dat',
#                   'run5_test/1D/SAXS/Reduction/sone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat')
# compare_dat_files('atT_smaller/OneD_1dd_WAXS_01/Reduction/Hor_all_WAXS.dat', 
#                   'run5_test/1D/WAXS/Reduction/b_tassone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat')

compare_dat_files('new_test/1D/SAXS/Reduction/sone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat',
                  'run5_test/1D/SAXS/Reduction/sone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat')
compare_dat_files('new_test/1D/WAXS/Reduction/b_tassone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat', 
                  'run5_test/1D/WAXS/Reduction/b_tassone_Hor_scan_Run5_RampT20_ctr0_scan1_0000.dat')

# What if the SAXS data comes in before the WAXS data