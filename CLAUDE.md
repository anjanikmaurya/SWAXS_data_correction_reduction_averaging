# SWAXS Data Processing Pipeline Documentation

## Overview

This codebase processes Small and Wide Angle X-ray Scattering (SAXS/WAXS) data from synchrotron experiments. The pipeline performs data reduction, normalization, and correction steps to convert raw detector images into scientifically meaningful scattering profiles.

## Project Structure

```
SWAXS_data_reduction_correction_Analysis/
├── Step1_SWAXS_Reduction_Normlisation_combined_macv4.ipynb  # Original notebook (legacy)
├── correction_integration_demo.py                          # NEW: Refactored data processing pipeline
├── config.yml                                              # Configuration file for refactored pipeline
├── Step2_SWAXS_bkg_subtract_01.ipynb                       # Background subtraction
├── Step3_SAXS_fit_porod_sub_v4.ipynb                       # SAXS analysis & fitting
├── Step3_WAXS_AmorphCal_and_Fitting_v3.ipynb              # WAXS crystallinity analysis
├── atT/                                                     # Experimental data directories (legacy format)
├── atT_smaller/                                             # Smaller experimental dataset 
├── run5_test/                                               # Test dataset for refactored implementation
├── python-version/                                          # Alternative Python implementation
└── pyproject.toml                                          # Project dependencies
```

## Data Correction Pipeline

### Refactored Pipeline (`correction_integration_demo.py`)

**Status: ACTIVE DEVELOPMENT - Core functionality implemented**

A significant refactoring of the original Step1 notebook has been implemented in `correction_integration_demo.py`. This modular approach provides the same data corrections with improved code organization and configuration management through an Experiment class.

**Key Components:**
- `Experiment` class (line 28-121): Configuration management and core functionality
- `Experiment.process_saxs_file()` (line 291-342): SAXS file processing with corrections and 1D integration
- `Experiment.process_waxs_file()` (line 344-395): WAXS file processing with corrections and 1D integration
- `calculate_correction_factors()` (line 219-263): Transmission and normalization factor calculation
- `read_raw_detector_file()` (line 155-183): Raw detector file reading
- `read_csv_parameters()` (line 122-153): CSV parameter extraction

**Key Improvements over Legacy:**
- **Multi-File Processing**: Processes all .raw files in directory structure
- **Experiment Class**: Centralized configuration management from `config.yml`
- **Modular Design**: Individual file processing methods for SAXS and WAXS
- **Error Handling**: Better validation and error messages
- **Flexible Output**: Maintains directory structure in output

### Legacy Pipeline (`Step1_SWAXS_Reduction_Normlisation_combined_macv4.ipynb`)

**Key Functions:**
- `process_SAXS_data()` (line 89-200): Processes raw SAXS detector files
- `process_WAXS_data()` (line 120-220): Processes raw WAXS detector files
- `calculate_sld_mu_thickness()` (line 45-75): Calculates material properties for normalization

**Critical Data Correction Steps (Both Implementations):**

1. **Detector Calibration & Geometry**
   - PyFAI integration using PONI files for SAXS (`atT_SAXS.poni`) and WAXS (`atT_WAXS.poni`)
   - Beam center coordinates: SAXS (491.215, 626.705), WAXS (668.647, 105.765)
   - Detector distances: SAXS 2.91m, WAXS 0.158m
   - **Refactored**: Parameters loaded from `config.yml` for flexibility

2. **Dark Current & Offset Correction**
   ```python
   # Legacy implementation (hardcoded)
   i0_corrected = i0_raw - i0_offset
   bstop_corrected = bstop_raw - bstop_offset
   
   # Refactored implementation (configurable)
   i0_corrected, bstop_corrected = correct_detector_offsets(
       i0_raw, bstop_raw, config['experimental']['i0_offset'], 
       config['experimental']['bstop_offset'])
   ```
   - Default offsets: `i0_offset = 0.285816667`, `bstop_offset = 0.477621111`

3. **Transmission & Thickness Normalization**
   ```python
   # Refactored approach with better error handling
   corrections = process_transmission_correction(
       i0_raw=parameters['i0_avg'],
       bstop_raw=parameters['bstop_avg'], 
       compound=config['compound'],
       config=config
   )
   ```
   - Uses configurable material properties (default: PE, density 0.92 g/cm³, energy 12 keV)
   - Enhanced Beer-Lambert law implementation with xraydb for absorption coefficients
   - Better error handling for invalid transmission values

4. **Mask Application**
   - SAXS: `RT_SAXS_mask_03.edf` (1043×981 pixels)
   - WAXS: `RT_WAXS_mask.edf` (195×487 pixels)
   - **Refactored**: Mask file paths configurable via `config.yml`

5. **1D Radial Integration**
   ```python
   # Both implementations use PyFAI integration
   q, I, error = ai.integrate1d(detector_image, npt_radial, 
                               error_model='poisson',
                               mask=mask, 
                               normalization_factor=corrections['normfactor'])
   ```
   - **Refactored**: Integration parameters configurable (npt_radial, error_model)

## Refactored Implementation Details

### Configuration Management System (`config.yml`)

The refactored implementation introduces a YAML-based configuration system managed by the Experiment class that addresses hardcoded parameter issues:

```yaml
# Example config.yml structure
data_directory: "run5_test/2D"  # Path to 2D data directory
poni_directory: "run5_test"     # Directory containing PONI and mask files  
compound: "C2H4"                # Material formula

detector_shapes:
  saxs: [1043, 981]
  waxs: [195, 487]

poni_files:
  saxs: "atT_SAXS.poni"
  waxs: "atT_WAXS.poni"

mask_files:
  saxs: "RT_SAXS_mask_03.edf"
  waxs: "RT_WAXS_mask.edf"

# Simplified experimental parameters (offsets set to 0 based on Step1 analysis)
energy_keV: 12
density_g_cm3: 0.92
i0_offset: 0
bstop_offset: 0
i0_air: 21.279793333
bstop_air: 18.65027
thickness: null  # Calculated from transmission if not provided

# Integration parameters
npt_radial: 1000
error_model: "poisson"
```

### File Structure Requirements

The refactored implementation expects a specific directory structure:

```
data_directory/
├── [experiment_dir]/
│   ├── SAXS/
│   │   ├── *.raw  # Raw detector files
│   │   └── *.csv  # Metadata files (same name as .raw)
│   └── WAXS/
│       ├── *.raw  # Raw detector files  
│       └── *.csv  # Metadata files (same name as .raw)
└── ...

Output: 1D/[experiment_dir]/
├── *.dat  # Processed SAXS/WAXS files (3-column: q, intensity, error)
└── ...
```

### Testing Requirements
**Priority: CRITICAL**

The refactored implementation requires comprehensive validation:

1. **Correctness Validation**
   - Compare q-vectors and intensities against original notebook results
   - Verify transmission factors and thickness calculations
   - Validate normalization factor computations
   - Check error propagation accuracy

2. **Configuration Testing**
   - Test with different material compounds
   - Verify detector parameter flexibility
   - Test with different experimental setups

3. **Error Handling**
   - Test file not found scenarios
   - Validate invalid transmission values
   - Check malformed configuration files

### Remaining Refactoring Requirements

### 1. **Multi-File Averaging**
**Priority: High - Currently Missing**

The refactored implementation processes individual files but lacks:
- Multi-file averaging capabilities from original notebook
- Time-series processing for multiple exposures  
- Statistical analysis of multiple measurements
- Error propagation in averaging calculations

### 2. **Flexible File Structure Support**
**Priority: Medium**

Future enhancements needed:
- Support for different file naming conventions
- Multiple data organization schemes  
- Automated metadata parsing from experimental files


## Key Dependencies

- `pyFAI`: Detector calibration and azimuthal integration
- `fabio`: Scientific image file I/O
- `periodictable`: Material property calculations
- `xraydb`: X-ray absorption coefficients
- `numpy`, `pandas`: Numerical processing
- `matplotlib`: Visualization

## Usage Notes

1. **Experimental Setup Requirements:**
   - PONI calibration files for each detector
   - Dark current measurements (shutter closed)
   - Empty cell/air measurements for background
   - Sample transmission measurements

2. **Data Quality Checks:**
   - Monitor I0 and beamstop stability over time
   - Verify transmission factors are physically reasonable
   - Check for detector artifacts in mask regions

3. **Processing Workflows:**

   **Refactored Pipeline (Multi-File Processing):**
   ```bash
   # Ensure config.yml is properly configured
   python correction_integration_demo.py
   
   # Then proceed with background subtraction and analysis
   jupyter notebook Step2_SWAXS_bkg_subtract_*.ipynb
   jupyter notebook Step3_*.ipynb
   ```
   
   **Legacy Pipeline (Multi-File with Averaging):**
   ```bash
   # Step 1: Raw data reduction with averaging
   jupyter notebook Step1_SWAXS_Reduction_*.ipynb
   
   # Step 2: Background subtraction
   jupyter notebook Step2_SWAXS_bkg_subtract_*.ipynb
   
   # Step 3: Analysis (SAXS fitting or WAXS crystallinity)
   jupyter notebook Step3_*.ipynb
   ```

## Important Notes for Refactored Implementation

**⚠️ CRITICAL WARNINGS:**

1. **Validation Required**: The refactored implementation (`correction_integration_demo.py`) is **UNDER ACTIVE DEVELOPMENT** and must be validated against known good results from the original notebook before production use.

2. **Configuration Dependency**: Requires a properly configured `config.yml` file with correct file paths and experimental parameters. Missing or incorrect configuration will cause processing failures.

3. **Different Data Structure**: The refactored implementation expects data in `2D/SAXS/` and `2D/WAXS/` subdirectories with CSV metadata files, different from the legacy `atT/` structure.

4. **Correction Factor Changes**: Based on Step1 notebook analysis, dark current offsets have been set to 0, which differs from earlier assumptions.

**Recommended Testing Protocol:**
1. Run both implementations on identical test data
2. Compare q-vectors, intensities, and error values numerically  
3. Verify transmission calculations and normalization factors
4. Test with different experimental conditions and materials
5. Validate error propagation and edge cases

This pipeline provides a robust foundation for SAXS/WAXS data processing. The original notebook remains the validated production system, while the refactored implementation offers improved modularity for future development.