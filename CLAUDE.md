# SWAXS Data Processing Pipeline Documentation

## Overview

This codebase processes Small and Wide Angle X-ray Scattering (SAXS/WAXS) data from synchrotron experiments. The pipeline performs data reduction, normalization, and correction steps to convert raw detector images into scientifically meaningful scattering profiles.

IMPORTANT: Don't test the code after running unless told otherwise
IMPORTANT: to run python code with imports, use uv run to use the virtual env

Also, you may find some diagnostic issues with the pt library on vscode, but these are not actual issues. You can ignore these.

## Project Structure

```
SWAXS_data_reduction_correction_Analysis/
├── src/                                                     # Main processing modules
│   ├── main_correction_reduction_v1.py                     # Core data processing pipeline
│   ├── process_metadata.py                                 # CSV/metadata processing utilities
│   ├── files_averaging.ipynb                              # Multi-file averaging functionality
│   └── utils/                                              # Utility modules
│       └── read_dat_metadata.py                           # .dat file parsing utilities
├── bin/                                                     # Command-line utilities and scripts
│   ├── compare_dat_files.py                               # Validation tool for comparing .dat outputs
│   └── copy_data_structure.py                             # Data organization utilities
├── config.yml                                              # Configuration file for processing pipeline
├── pyproject.toml                                          # Project dependencies
├── run5_test/                                               # Test dataset for validation
└── [legacy notebooks will be moved/archived]               # Step1-3 notebooks (transitioning out)
```

## Data Correction Pipeline

### Core Processing Module (`src/main_correction_reduction_v1.py`)

**Status: ACTIVE DEVELOPMENT - Core functionality implemented**

The main data processing pipeline has been refactored into a modular Python implementation in `src/`. This provides the same data corrections as the original notebooks with improved code organization and configuration management.

**Key Components:**
- `Experiment` class: Configuration management and core functionality
- `Experiment.process_saxs_file()`: SAXS file processing with corrections and 1D integration
- `Experiment.process_waxs_file()`: WAXS file processing with corrections and 1D integration
- `calculate_correction_factors()`: Transmission and normalization factor calculation
- `read_raw_detector_file()`: Raw detector file reading

**Supporting Modules:**
- `src/process_metadata.py`: CSV/metadata processing utilities
  - `process_csv_metadata()`: Extracts i0, bstop values from experimental CSV files
  - `find_row_number_to_read()`: Maps raw file indices to CSV metadata rows
- `src/utils/read_dat_metadata.py`: .dat file parsing utilities
  - `read_dat_data_metadata()`: Parses processed .dat files with embedded metadata
- `src/files_averaging.ipynb`: Multi-file averaging functionality (work in progress)

**Key Improvements over Legacy:**
- **Modular Architecture**: Separated concerns across multiple focused modules
- **Configuration Management**: Centralized YAML-based configuration in `config.yml`
- **Multi-File Processing**: Processes all .raw files in directory structure automatically
- **Better Error Handling**: Comprehensive validation and error messages
- **Flexible Output**: Maintains directory structure with embedded metadata

### Command-Line Utilities (`bin/`)

**Validation and Data Management Tools:**
- `bin/compare_dat_files.py`: Numerical comparison tool for validating .dat file outputs
  - `compare_dat_files()`: Compares processed data files within specified tolerance
  - `write_comparison_results()`: Generates detailed comparison reports
- `bin/copy_data_structure.py`: Data organization and structure management utilities

## Critical Data Correction Steps

The pipeline implements essential SWAXS data corrections:

1. **Detector Calibration & Geometry**
   - PyFAI integration using PONI files for SAXS and WAXS detectors
   - Configurable beam center coordinates and detector distances
   - Parameters loaded from `config.yml` for experimental flexibility

2. **Dark Current & Offset Correction**
   ```python
   # Modular implementation with configurable offsets
   i0_corrected, bstop_corrected = correct_detector_offsets(
       i0_raw, bstop_raw,
       config['experimental']['i0_offset'],
       config['experimental']['bstop_offset'])
   ```

3. **Transmission & Thickness Normalization**
   ```python
   # Enhanced approach with material property calculations
   corrections = process_transmission_correction(
       i0_raw=parameters['i0_avg'],
       bstop_raw=parameters['bstop_avg'],
       compound=config['compound'],
       config=config
   )
   ```
   - Configurable material properties (density, energy, composition)
   - Enhanced Beer-Lambert law implementation using xraydb
   - Comprehensive error handling for invalid transmission values

4. **Mask Application**
   - Detector-specific mask files for artifact removal
   - Configurable mask file paths via `config.yml`
   - Different masks for SAXS and WAXS detectors

5. **1D Radial Integration**
   ```python
   # PyFAI-based azimuthal integration
   q, I, error = ai.integrate1d(detector_image, npt_radial,
                               error_model='poisson',
                               mask=mask,
                               normalization_factor=corrections['normfactor'])
   ```
   - Configurable integration parameters (npt_radial, error_model)
   - Poisson error model for proper uncertainty propagation

## Implementation Details

### Configuration Management System (`config.yml`)

The modular implementation uses a YAML-based configuration system managed by the Experiment class:

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

The modular implementation expects a specific directory structure:

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

### Testing and Validation

**Priority: CRITICAL**

The modular implementation requires comprehensive validation:

1. **Correctness Validation**
   - Compare q-vectors and intensities using `bin/compare_dat_files.py`
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

**Validation Tools:**
- `bin/compare_dat_files.py`: Numerical comparison of processed .dat files
- `src/utils/read_dat_metadata.py`: Metadata extraction for validation

### Development Roadmap

**Current Priorities:**

1. **Multi-File Averaging** (`src/files_averaging.ipynb`)
   - **Status**: Work in Progress
   - **Needs**: Complete integration with main processing pipeline
   - **Features**: Time-series processing, statistical analysis, error propagation

2. **Enhanced Validation Tools**
   - **Status**: Partial implementation in `bin/`
   - **Needs**: Comprehensive test suite integration
   - **Features**: Automated regression testing, numerical validation

3. **Flexible File Structure Support**
   - **Status**: Future enhancement
   - **Features**: Multiple naming conventions, automated metadata parsing


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

   **Main Processing Pipeline:**
   ```bash
   # Ensure config.yml is properly configured
   uv run src/main_correction_reduction_v1.py

   # Validate results (optional)
   uv run bin/compare_dat_files.py file1.dat file2.dat

   # Then proceed with background subtraction and analysis
   jupyter notebook Step2_SWAXS_bkg_subtract_*.ipynb
   jupyter notebook Step3_*.ipynb
   ```

   **Development and Testing:**
   ```bash
   # Data validation
   uv run bin/compare_dat_files.py reference.dat output.dat

   # Data structure management
   uv run bin/copy_data_structure.py

   # Multi-file averaging (when completed)
   jupyter notebook src/files_averaging.ipynb
   ```

## Important Notes

**⚠️ CRITICAL WARNINGS:**

1. **Validation Required**: The modular implementation (`src/main_correction_reduction_v1.py`) is **UNDER ACTIVE DEVELOPMENT** and must be validated against known good results before production use.

2. **Configuration Dependency**: Requires a properly configured `config.yml` file with correct file paths and experimental parameters. Missing or incorrect configuration will cause processing failures.

3. **Data Structure Requirements**: The implementation expects data in `2D/SAXS/` and `2D/WAXS/` subdirectories with corresponding CSV metadata files.

4. **Correction Factor Updates**: Dark current offsets have been revised based on systematic analysis.

**Recommended Testing Protocol:**
1. Use `bin/compare_dat_files.py` for numerical validation
2. Compare q-vectors, intensities, and error values against reference data
3. Verify transmission calculations and normalization factors
4. Test with different experimental conditions and materials
5. Validate error propagation and edge cases

**Development Focus:**
- The `src/` directory contains the core processing logic
- The `bin/` directory provides validation and utility tools
- Legacy notebooks are being phased out in favor of modular Python implementation

This modular pipeline provides a robust foundation for SWAXS data processing with improved maintainability and extensibility.