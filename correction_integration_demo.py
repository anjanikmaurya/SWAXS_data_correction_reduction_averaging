#!/usr/bin/env python3
"""
SWAXS Data Correction and 1D Integration Pipeline
=================================================
Refactored implementation of Step1 notebook with improved modularity and configuration management.
Processes all .raw files in directory structure with essential data corrections and 1D integration.

Status: ACTIVE DEVELOPMENT - Core functionality implemented with accurate calculations, requires extension

Key Features:
- Experiment class for centralized configuration management from config.yml
- Multi-file processing (all .raw files in SAXS/WAXS subdirectories)
- Essential data corrections (transmission, normalization, material properties)
- Accurate SLD calculation using periodictable for precise material properties
- 1D radial integration with PyFAI
- Configurable experimental parameters and detector geometry
- Output in .dat format compatible with existing analysis workflows
- Uses PDI files instead of CSV

File Structure Requirements:
- Input: {data_directory}/[experiment_dirs]/SAXS/*.raw and WAXS/*.raw
- Metadata: Corresponding .csv files with same name as .raw files
- Output: 1D/[experiment_dirs]/*.dat (3-column: q, intensity, error)

Logging rules:
1. Logs all config information just once — when Experiment._load_and_assign_config() is called
2. Also performs logging for all corrections called on the very first run of SAXS and the first run of WAXS (if both are run).

I should first complete the logic and then afterwards shard my directories.
When completing the other logic, I should probably add it to another file
"""
# TODO: Verify that it is working with the current data
# TODO: 1-5: Testing with different 1-5 data.
# TODO: 17-2: Look at how to add these into
# TODO: integrate with different steps
# TODO: Switch to pathlib

# TODO: Split up into different files 
# TODO: Logging -- add PONI files, function by function
# TODO: Plots: Use Anjani's function to create y offsets
# TODO: Read CSVs. If multiple rows in CSV, then index of CSV corresponds to number after scan
# TODO: Specify
import os
import glob
import numpy as np
import yaml
import fabio
import pyFAI
import periodictable as pt
import xraydb
from typing import Dict, List, Tuple
import re
from pathlib import Path
import logging
import pandas as pd

# Importing files
import utils
import process_pdi_metadata
import process_csv_metadata

# Start logger, overwrites existing logging files, includes timestamp in message
logging.basicConfig(filename = "corrections.log",
                    filemode = "w", 
                    level = logging.INFO,
                    format="%(asctime)s -  %(name)s - %(levelname)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger('swaxs_pipeline')
class Experiment:
    # Type hints for dynamically assigned config attributes
    data_directory: str
    poni_directory: Path
    compound: str
    energy_keV: float
    density_g_cm3: float
    detector_shapes: Dict[str, List[int]]
    poni_files: Dict[str, str]
    mask_files: Dict[str, str]
    read_metadata_function: str
    thickness: float
    npt_radial: int
    error_model: str
    i0_offset: float
    i0_air: float
    bstop_offset: float
    bstop_air: float
    mode: str
    metadata_format: str
    """
    SWAXS Experiment class for managing configuration and data processing.
    
    Loads experimental parameters from config.yml and provides methods for
    processing SAXS/WAXS data with proper corrections and normalization.
    """
    
    def __init__(self, config_path: str = "config.yml"):
        """
        Initialize Experiment with configuration from YAML file.
        
        Parameters
        ----------
        config_path : str
            Path to the configuration file
        """
        self.config_path = config_path
        # Logging state tracking for first-run detailed logging
        self._saxs_log_pending = True
        self._waxs_log_pending = True
        
        self._load_and_assign_config()
        self._setup_directories()
        self._load_pyfai_integrators()
        
    def _load_and_assign_config(self):
        """Load configuration from YAML file and assign all values to instance variables."""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Log all configuration information once
        logger.info("_load_and_assign_config called")

        logger.info("=" * 60)
        logger.info("EXPERIMENT CONFIGURATION")
        logger.info("=" * 60)
        logger.info(f"Configuration file: {self.config_path}")
        logger.info("")
        
        # Loop through all config keys and assign them as instance variables
        for key, value in config.items():
            setattr(self, key, value)
            # Log each configuration parameter
            if isinstance(value, dict):
                logger.info(f"{key}:")
                for sub_key, sub_value in value.items():
                    logger.info(f"  {sub_key}: {sub_value}")
            elif isinstance(value, list):
                logger.info(f"{key}: {value}")
            else:
                logger.info(f"{key}: {value}")
        
        # Make mode upper-case (SAXS, SWAXS, WAXS) for case-insensitive string comparison
        self.mode = self.mode.upper()
        logger.info("")
        logger.info("=" * 60)
    
    def _setup_directories(self):
        """Setup computed directory paths."""
        self.data_directory_2d = Path(self.data_directory) / "2D"
        self.poni_directory = Path(self.data_directory) / Path(self.poni_directory)
        self.output_directory_1d = Path(self.data_directory) / "1D"
        self.output_directory_1d.mkdir(parents=True, exist_ok=True)
        
        self.saxs_subdir = "SAXS" 
        self.waxs_subdir = "WAXS"
    
    def _load_saxs_integrator(self):
        """Load PyFAI integrator and mask for SAXS detector."""
        self.saxs_poni_path = self.poni_directory / self.poni_files['saxs']
        
        if not self.saxs_poni_path.exists():
            raise FileNotFoundError(f"SAXS PONI file not found: {self.saxs_poni_path}")
            
        self.ai_saxs = pyFAI.load(str(self.saxs_poni_path))
        logger.info(f"Loaded SAXS PyFAI integrator from: {self.saxs_poni_path}")
        
        # Load SAXS detector mask
        if self.mask_files['saxs'] is not None:
            saxs_mask_path = self.poni_directory / self.mask_files['saxs']
            self.saxs_mask = fabio.open(str(saxs_mask_path)).data
            logger.info(f"Loaded SAXS mask from: {saxs_mask_path}")
        else:
            self.saxs_mask = None
            logger.info("No SAXS mask file specified")
    
    def _load_waxs_integrator(self):
        """Load PyFAI integrator and mask for WAXS detector."""
        self.waxs_poni_path = self.poni_directory / self.poni_files['waxs']
        
        if not self.waxs_poni_path.exists():
            raise FileNotFoundError(f"WAXS PONI file not found: {self.waxs_poni_path}")
            
        self.ai_waxs = pyFAI.load(str(self.waxs_poni_path))
        logger.info(f"Loaded WAXS PyFAI integrator from: {self.waxs_poni_path}")
        
        # Load WAXS detector mask
        if self.mask_files['waxs'] is not None:
            waxs_mask_path = self.poni_directory / self.mask_files['waxs']
            self.waxs_mask = fabio.open(str(waxs_mask_path)).data 
            logger.info(f"Loaded WAXS mask from: {waxs_mask_path}")
        else:
            self.waxs_mask = None
            logger.info("No WAXS mask file specified")

    def _load_pyfai_integrators(self):
        """Load PyFAI integrator objects based on processing mode."""
        
        if self.mode not in ['SAXS', 'WAXS', 'SWAXS']:
            raise ValueError(f"Invalid mode '{self.mode}'. Must be 'SAXS', 'WAXS', or 'SWAXS'")

        if self.mode == 'WAXS':
            self._load_waxs_integrator()
            self._saxs_log_pending = False
            
        elif self.mode == 'SAXS':
            self._load_saxs_integrator()
            self._waxs_log_pending = False

        else:
            self._load_saxs_integrator()
            self._load_waxs_integrator()


    def read_raw_detector_15_image(self, raw_file_path: Path, detector_type: str) -> np.ndarray:
        """
        Read raw detector file and return detector image array.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the .raw detector file
        detector_type : str
            Type of detector ('SAXS' or 'WAXS')
            
        Returns
        -------
        np.ndarray
            Detector image array with proper shape
        """
        # Log function call for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"Called read_raw_detector_15_image() for {detector_type} detector")
        
        if detector_type.upper() == 'SAXS':
            shape = tuple(self.detector_shapes['saxs'])
        elif detector_type.upper() == 'WAXS':
            shape = tuple(self.detector_shapes['waxs'])
        else:
            raise ValueError("detector_type must be 'SAXS' or 'WAXS'")
        
        if not raw_file_path.exists():
            raise FileNotFoundError(f"Raw detector file not found: {raw_file_path}")
        
        data = np.fromfile(str(raw_file_path), dtype=np.int32).reshape(shape)
        
        # Log data information for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"  Raw data loaded: shape {data.shape}, dtype {data.dtype}")
            logger.info(f"  Data range: min={data.min()}, max={data.max()}, mean={data.mean():.2f}")
        
        return data

    def calculate_sld_mu_thickness(self, transmission: float) -> Dict[str, float]:
        """
        Calculate material properties from transmission using experiment parameters.
        
        Parameters
        ----------
        transmission : float
            Measured transmission value
            
        Returns
        -------
        Dict[str, float]
            Dictionary containing mu, sld, and thickness_m
        """
        # Log function call for first runs only
        # if ((self._saxs_log_pending) or (self._waxs_log_pending)):
        #     logger.info(f"Called calculate_sld_mu_thickness() with transmission = {transmission:.6f}")
        
        if transmission <= 0 or transmission > 1:
            raise ValueError("Transmission must be in the range (0, 1]")
        
        # Use xraydb for absorption coefficient calculation
        mu_xraydb = xraydb.material_mu(self.compound, 
                                     energy=self.energy_keV * 1000, 
                                     density=self.density_g_cm3) * 100  # Convert to 1/m
        
        # Calculate scattering length density (SLD) using periodictable
        material = pt.formula(self.compound) # This works properly, in spite of IDE doubting it
        sld, mu_pt = pt.xray_sld(material, energy=self.energy_keV, density=self.density_g_cm3)

        
        # Calculate thickness: T = exp(-mu * t) → t = -ln(T) / mu
        thickness = -np.log(transmission) / mu_xraydb
        
        return {
            "mu": mu_xraydb, 
            "sld": sld, 
            "thickness_m": thickness
        }

    def create_output_directory(self, raw_file_path: Path, detector_type: str) -> Path:
        """
        TODO: Move this to the other creation directory
        Create output directory structure with SAXS/Reduction or WAXS/Reduction subdirectories.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the raw file being processed
        detector_type : str
            Type of detector ('SAXS' or 'WAXS')
            
        Returns
        -------
        str
            Output directory path
        """
        # Log function call for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"Called create_output_directory() for {detector_type} detector")
            
        output_dir = self.output_directory_1d / detector_type.upper() / "Reduction"
        output_dir.mkdir(parents=True, exist_ok=True)
    
        # Log output directory for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"  Output directory: {output_dir}")
            
        return output_dir


    def get_corrections_full(self, raw_file_path: Path, detector_type: str):
        """Processed metadata file and computes correction factors"""
        # Log function call for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"Called get_corrections_full() for {detector_type} detector")
        
        # Use specified metadata function from utils if available, otherwise use PDI processing
        if self.read_metadata_function is not None and hasattr(utils, self.read_metadata_function):
            metadata_function = getattr(utils, self.read_metadata_function)
            i0, bstop, metadata_path = metadata_function(raw_file_path)
        else:
            # Using current PDI metadata processing
            if self.metadata_format == "csv":
                i0, bstop, metadata_path = process_csv_metadata.process_csv_metadata(raw_file_path)
            elif self.metadata_format == "pdi":
                i0, bstop, metadata_path = process_pdi_metadata.process_pdi_full(raw_file_path, detector_type)
            else:
                raise RuntimeError(f"Metadata format {self.metadata_format} not support. Needs to be csv or pdi.")
        
        # Log metadata processing for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"  Processing metadata from: {metadata_path}")
            logger.info(f"  Raw I0: {i0:.6f}, Raw Bstop: {bstop:.6f}")
        
        # Calculate correction factors
        transmission_factor_raw = bstop
        transmission_ratio = bstop / i0
        
        # Log transmission calculations for first runs only
        if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
            (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
            logger.info(f"  Calculating transmission: ratio = {transmission_ratio:.6f}")
        
        # Calculate material properties if thickness not provided
        if self.thickness is None:
            
            material_props = self.calculate_sld_mu_thickness(transmission_ratio)
            thickness = material_props['thickness_m']
            
            # Log material calculations for first runs only
            if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
                (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
                logger.info(f"  Material properties calculated:")
                logger.info(f"    Compound: {self.compound}")
                logger.info(f"    Energy: {self.energy_keV} keV")
                logger.info(f"    Density: {self.density_g_cm3} g/cm³")
                logger.info(f"    Absorption coefficient: {material_props['mu']:.6f} 1/m")
                logger.info(f"    SLD: {material_props['sld']:.6e}")
                logger.info(f"    Calculated thickness: {thickness:.6f} m")
        else:
            thickness = self.thickness
            if ((detector_type.upper() == 'SAXS' and self._saxs_log_pending) or 
                (detector_type.upper() == 'WAXS' and self._waxs_log_pending)):
                logger.info(f"  Using provided thickness: {thickness:.6f} m")
        
        # Calculate normalization factor: trans_factor * i0_corrected
        # Based on Step1 notebook: normfactor = trans_factor*i0_avg_all[0]
        normalization_factor = transmission_factor_raw * i0
        
        return {
            'i0_corrected': i0,
            'bstop_corrected': bstop,
            'transmission_factor': transmission_factor_raw,
            'transmission_ratio': transmission_ratio,
            'thickness': thickness,
            'normalization_factor': normalization_factor,
            'metadata_path': metadata_path
        }
        
    def process_saxs_file(self, raw_file_path: Path):
        """
        Process a single SAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the SAXS .raw file
            
        Returns
        -------
        Dict
            Dictionary containing plotting data: q, intensity, error, filename, corrections
        """
                
        # Log detailed information for first SAXS file only
        if self._saxs_log_pending:
            logger.info("")
            logger.info("=" * 60)
            logger.info("FIRST SAXS FILE PROCESSING")
            logger.info("=" * 60)
            logger.info(f"Raw file: {raw_file_path}")
            logger.info(f"SAXS PONI file: {self.saxs_poni_path}")
            logger.info(f"SAXS mask file: {self.poni_directory / self.mask_files['saxs'] if self.mask_files['saxs'] else 'None'}")
            logger.info("")
        
        # Read raw detector data
        detector_data = self.read_raw_detector_15_image(raw_file_path, 'SAXS')
        
        corrections = self.get_corrections_full(raw_file_path, "SAXS")
        
        # Log corrections for first SAXS file only
        if self._saxs_log_pending:
            logger.info(f"  Metadata file: {corrections['metadata_path']}")

            self._saxs_log_pending = False
        
        # Create output directory and filename
        # TODO: Make output directory only once.
        output_dir = self.create_output_directory(raw_file_path, 'SAXS')
        output_filename = raw_file_path.name.replace('.raw', '.dat')
        output_filename = output_filename.removeprefix("sone_")

        output_path = output_dir / output_filename
        output_path = output_path.with_stem(f"{output_path.stem}_SAXS")
        # Perform 1D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_saxs.integrate1d(
            detector_data, 
            self.npt_radial,
            error_model=self.error_model,
            mask=self.saxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=str(output_path)
        )
                
        add_metadata_to_dat(dat_file_path= output_path, metadata_file_path=corrections['metadata_path'])
        # Return plotting data
        return {
            'q': q,
            'intensity': intensity,
            'error': error,
            'filename': output_filename,
            'corrections': corrections,
            'raw_file_path': raw_file_path
        }
    def process_waxs_file(self, raw_file_path: Path):
        """
        Process a single WAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the WAXS .raw file
            
        Returns
        -------
        Dict
            Dictionary containing plotting data: q, intensity, error, filename, corrections
        """        
        
        # Log detailed information for first WAXS file only
        if self._waxs_log_pending:
            logger.info("")
            logger.info("=" * 60)
            logger.info("FIRST WAXS FILE PROCESSING")
            logger.info("=" * 60)
            logger.info(f"Raw file: {raw_file_path}")
            logger.info(f"WAXS PONI file: {self.waxs_poni_path}")
            logger.info(f"WAXS mask file: {self.poni_directory / self.mask_files['waxs'] if self.mask_files['waxs'] else 'None'}")
            logger.info(f"Detector shape: {self.detector_shapes['waxs']}")
            logger.info("")
        
        detector_data = self.read_raw_detector_15_image(raw_file_path, 'WAXS')
        corrections = self.get_corrections_full(raw_file_path, "WAXS")

        # Log corrections for first WAXS file only
        if self._waxs_log_pending:
            logger.info(f"  Metadata file: {corrections['metadata_path']}")
            logger.info("=" * 60)
            self._waxs_log_pending = False

        # Create output directory and filename
        output_filename = raw_file_path.name.replace('.raw', '.dat')
        output_filename = output_filename.removeprefix("b_tassone_")
        
        output_dir = self.create_output_directory(raw_file_path, 'WAXS')
        output_path = output_dir / output_filename
        output_path = output_path.with_stem(f"{output_path.stem}_WAXS")
        
        # Perform D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_waxs.integrate1d(
            detector_data,
            self.npt_radial,
            error_model=self.error_model,
            mask=self.waxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=str(output_path)
        )
        add_metadata_to_dat(dat_file_path= output_path, metadata_file_path=corrections['metadata_path'])

        # Return plotting data
        return {
            'q': q,
            'intensity': intensity,
            'error': error,
            'filename': output_filename,
            'corrections': corrections,
            'raw_file_path': raw_file_path
        }


def find_all_raw_files(experiment: Experiment, data_directory_path: Path) -> Tuple[List[Path], List[Path]]:
    
    """
    Find all .raw files in SAXS and WAXS directories.
    
    Parameters
    ----------
    data_directory_path : Path
        Full path to the data directory (e.g., "run12_test/2D/")
        
    Returns
    -------
    Tuple[List[Path], List[Path]]
        Lists of SAXS and WAXS .raw file paths
    """
    if (experiment.mode != "WAXS"):
        saxs_pattern = str(data_directory_path / "**/SAXS/**/*.raw")
        saxs_files = [Path(f) for f in glob.glob(saxs_pattern, recursive=True) 
                    if not f.endswith('.raw.pdi')]
    else:
        # Don't include SAXS files in just WAXS
        saxs_files = []
    
    if (experiment.mode != "SAXS"):
        waxs_pattern = str(data_directory_path / "**/WAXS/**/*.raw")
        waxs_files = [Path(f) for f in glob.glob(waxs_pattern, recursive=True)
                    if not f.endswith('.raw.pdi')]
    else:
        # Don't include WAXS files in just SAXS
        waxs_files = []
        
    return sorted(saxs_files), sorted(waxs_files)

def full_correction_integration(plotting=False):
    """
    Main function to process all .raw files in the dataset.
    
    Parameters
    ----------
    plotting : bool, default False
        If True, return plotting data instead of just processing files
        
    Returns
    -------
    Dict or None
        If plotting=True, returns dictionary with plotting data.
        If plotting=False, returns None (just processes files).
    """
    print("SWAXS Data Correction and 1D Integration Demo")
    print("=" * 50)
    
    # Initialize experiment from config
    experiment = Experiment("config.yml")
    print(f"Loaded configuration for compound: {experiment.compound}")
    print(f"Energy: {experiment.energy_keV} keV, Density: {experiment.density_g_cm3} g/cm³")
    print()
    
    # Find all .raw files to process
    saxs_files, waxs_files = find_all_raw_files(experiment, experiment.data_directory_2d)
    
    print(f"Found {len(saxs_files)} SAXS files and {len(waxs_files)} WAXS files to process")
    print()
    
    # Initialize plotting data collection if requested
    if plotting:
        plotting_data = {
            'saxs_data': [],
            'waxs_data': [],
            'metadata': {
                'timers': [],
                'temperatures': [], 
                'i0_values': [],
                'bstop_values': [],
                'filenames': []
            }
        }
    
    # Process all SAXS files
    if saxs_files:
        print("Processing SAXS files...")
        for saxs_file in saxs_files:
            saxs_result = experiment.process_saxs_file(saxs_file)
            if plotting:
                plotting_data['saxs_data'].append(saxs_result)
                # Extract metadata for monitoring plots
                corrections = saxs_result['corrections']
                plotting_data['metadata']['i0_values'].append(corrections['i0_corrected'])
                plotting_data['metadata']['bstop_values'].append(corrections['bstop_corrected'])
                plotting_data['metadata']['filenames'].append(saxs_result['filename'])
        print()
    
    # Process all WAXS files
    if waxs_files:
        print("Processing WAXS files...")
        for waxs_file in waxs_files:
            waxs_result = experiment.process_waxs_file(waxs_file)
            if plotting:
                plotting_data['waxs_data'].append(waxs_result)
        print()
    
    print("Processing complete!")
    
    # Return plotting data if requested
    if plotting:
        return plotting_data
    else:
        return None


def add_metadata_to_dat(dat_file_path: Path, metadata_file_path: str):
    """
    Add PDI metadata to the end of a .dat file as commented information.
    
    Parameters
    ----------
    dat_file_path : str
        Path to the .dat file to modify
    pdi_file_path : str
        Path to the .pdi file containing metadata
        
    Notes
    -----
    - Appends PDI data as comments (lines starting with #) to end of .dat file
    - If a line already starts with #, adds another # prefix  
    - Preserves original .dat file structure and data columns
    """
    if not dat_file_path.exists():
        raise FileNotFoundError(f"DAT file not found: {dat_file_path}")
    
    if not Path(metadata_file_path).exists():
        raise FileNotFoundError(f"PDI file not found: {metadata_file_path}")
    
    # Read the PDI file
    with open(metadata_file_path, 'r') as f:
        pdi_lines = f.readlines()
    
    # Append PDI data to the .dat file
    with open(dat_file_path, 'a') as f:
        f.write(f"\n#  Metadata from: {Path(metadata_file_path).name}\n")
        
        for line in pdi_lines:
            clean_line = line.rstrip()
            if clean_line:  # Skip empty lines
                # Add # prefix, and add another # if line already starts with #
                if clean_line.startswith('#'):
                    f.write(f"#{clean_line}\n")
                else:
                    f.write(f"# {clean_line}\n")
        
def main():
    return full_correction_integration()

if __name__ == "__main__":
    main()
