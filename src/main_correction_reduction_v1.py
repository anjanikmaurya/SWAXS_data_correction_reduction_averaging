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
# TODO: 17-2: Look at how to add these into
# TODO: Move reductions.log to poni directory
import glob
import numpy as np
import yaml
import fabio
import pyFAI
import periodictable as pt
import xraydb
from typing import Dict, List, Tuple
from pathlib import Path
import logging
import argparse

import process_metadata


# Logging will be configured after loading config in Experiment.__init__
logger = logging.getLogger('swaxs_pipeline')
# Suppress PyFAI warnings to reduce log spam while keeping error messages
logging.getLogger('pyFAI').setLevel(logging.ERROR)

class Experiment:
    # Type hints for dynamically assigned config attributes.
    data_directory: Path
    poni_directory: Path # Name of poni directory
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
    
    poni_path: Path # Path to poni directory
    _log_pending: bool # Specifies if the first-run log is currently pending

    def __init__(self, config_path: str):
        """
        Initialize Experiment with configuration from YAML file.
        
        Parameters
        ----------
        config_path : str
            Path to the configuration file
        """
        self.config_path = config_path
        self._log_pending = True

        self._load_and_assign_config()
        self._setup_directories()
        self._load_pyfai_integrators()

    def _load_and_assign_config(self):
        """Load configuration from YAML file and assign all values to instance variables."""
        
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        self.data_directory = Path(config['data_directory'])
        del config['data_directory'] # Don't reread the data directory -- important to not store as a string
        if not self.data_directory.exists():
            raise RuntimeError("Data Directory Not Found!")
        self._setup_logging()
        # Log all configuration information onc
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
    def _setup_logging(self):
        logging_path = self.data_directory / "reductions.log"
        logging.basicConfig(filename = logging_path,
            filemode = "w", 
            level = logging.INFO,
            format="%(asctime)s -  %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S")
        logger.info("_load_and_assign_config called")

        logger.info("=" * 60)
        logger.info("EXPERIMENT CONFIGURATION")
        logger.info("=" * 60)
    
    def _setup_directories(self):
        """Setup computed directory paths."""
        self.data_directory_2d = self.data_directory / "2D"
        if not self.data_directory_2d.exists():
            raise RuntimeError(f"No 2D Directory Found in {self.data_directory}. Please move your data to a 2D directory.")
        self.poni_path = self.data_directory / self.poni_directory
        if not self.poni_directory:
            raise RuntimeError(f"PONI Directory not found: ")

        self.output_directory_1d = Path(self.data_directory) / "1D"
        self.output_directory_1d.mkdir(parents=True, exist_ok=True)
        
        self.saxs_subdir = "SAXS" 
        self.waxs_subdir = "WAXS"
    
    def _load_saxs_integrator(self):
        """Load PyFAI integrator and mask for SAXS detector."""
        self.saxs_poni_path = self.poni_path / self.poni_files['saxs']
        
        if not self.saxs_poni_path.exists():
            raise FileNotFoundError(f"SAXS PONI file not found: {self.saxs_poni_path}")
            
        self.ai_saxs = pyFAI.load(str(self.saxs_poni_path))
        logger.info(f"Loaded SAXS PyFAI integrator from: {self.saxs_poni_path}")
        
        # Load SAXS detector mask
        if self.mask_files['saxs'] is not None:
            saxs_mask_path = self.poni_path / self.mask_files['saxs']
            self.saxs_mask = fabio.open(str(saxs_mask_path)).data
            logger.info(f"Loaded SAXS mask from: {saxs_mask_path}")
        else:
            self.saxs_mask = None
            logger.info("No SAXS mask file specified")
    
    def _load_waxs_integrator(self):
        """Load PyFAI integrator and mask for WAXS detector."""
        self.waxs_poni_path = self.poni_path / self.poni_files['waxs']
        
        if not self.waxs_poni_path.exists():
            raise FileNotFoundError(f"WAXS PONI file not found: {self.waxs_poni_path}")
            
        self.ai_waxs = pyFAI.load(str(self.waxs_poni_path))
        logger.info(f"Loaded WAXS PyFAI integrator from: {self.waxs_poni_path}")
        
        # Load WAXS detector mask
        if self.mask_files['waxs'] is not None:
            waxs_mask_path = self.poni_path / self.mask_files['waxs']
            self.waxs_mask = fabio.open(str(waxs_mask_path)).data 
            logger.info(f"Loaded WAXS mask from: {waxs_mask_path}")
        else:
            self.waxs_mask = None
            logger.info("No WAXS mask file specified")

    def _load_pyfai_integrators(self):
        """Load PyFAI integrator objects based on processing mode."""
        
        if self.mode not in ['SAXS', 'WAXS', 'SWAXS']:
            raise ValueError(f"Invalid mode '{self.mode}' in configuration file. Must be 'SAXS', 'WAXS', or 'SWAXS'")

        if self.mode == 'WAXS':
            self._load_waxs_integrator()
            
        elif self.mode == 'SAXS':
            self._load_saxs_integrator()

        else:
            self._load_saxs_integrator()
            self._load_waxs_integrator()


    def read_raw_detector_15_image(self, raw_file_path: Path, detector_type: str) -> np.ndarray:
        """
        Read raw detector file and return detector image array.
        
        Parameters
        ----------
        raw_file_path Path to the .raw detector file
        detector_type : Type of detector ('SAXS' or 'WAXS')
            
        Returns
        -------
        np.ndarray
            Detector image array with proper shape
        """
        # Log function call for first runs only
        if self._log_pending:
            logger.info(f"Called read_raw_detector_15_image() for {detector_type} detector")
        
        if detector_type.upper() == 'SAXS':
            shape = tuple(self.detector_shapes['saxs'])
        elif detector_type.upper() == 'WAXS':
            shape = tuple(self.detector_shapes['waxs'])

        data = np.fromfile(str(raw_file_path), dtype=np.int32).reshape(shape)

        
        if not raw_file_path.exists():
            raise FileNotFoundError(f"Raw detector file not found: {raw_file_path}")
        
        data = np.fromfile(str(raw_file_path), dtype=np.int32).reshape(shape)
        
        # Log data information for first runs only
        if self._log_pending:
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
        if self._log_pending:
            logger.info(f"Called calculate_sld_mu_thickness() with transmission = {transmission:.6f}")
        
        if transmission <= 0 or transmission > 1:
            raise ValueError("Transmission must be in the range (0, 1]")
        
        # Use xraydb for absorption coefficient calculation
        mu_xraydb = xraydb.material_mu(self.compound, 
                                     energy=self.energy_keV * 1000, 
                                     density=self.density_g_cm3) * 100  # Convert to 1/m
        
        # Calculate scattering length density (SLD) using periodictable
        material = pt.formula(self.compound) # pt library works properly even with red squiggle
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
        if self._log_pending:
            logger.info(f"Called create_output_directory() for {detector_type} detector")
            
        output_dir = self.output_directory_1d / detector_type.upper() / "Reduction"
        output_dir.mkdir(parents=True, exist_ok=True)
    
        # Log output directory for first runs only
        if self._log_pending:
            logger.info(f"  Output directory: {output_dir}")
            
        return output_dir


    def get_corrections_full(self, raw_file_path: Path, detector_type: str):
        """Processed metadata file and computes correction factors"""
        # Log function call for first runs only
        if self._log_pending:
            logger.info(f"Called get_corrections_full() for {detector_type} detector with file path: {raw_file_path}")
        
        # Use specified metadata function from utils if available, otherwise use CSV processing
            # Using current PDI metadata processing
        if self.metadata_format == "csv":
            metadata_dict = process_metadata.process_csv_metadata(raw_file_path)
        elif self.metadata_format == "pdi":
            # TODO
            metadata_dict = process_metadata.process_csv_metadata(raw_file_path)

        else:
            raise RuntimeError(f"Metadata format {self.metadata_format} not supported. Needs to be csv.")
        
        i0 = metadata_dict['i0']
        bstop = metadata_dict['bstop']
        if self._log_pending:
            logger.info(f"  Raw I0: {i0:.6f}, Raw Bstop: {bstop:.6f}")
        
        # Calculate correction factors
        transmission_factor_raw = bstop
        transmission_ratio = bstop / i0
        i0_corrected = i0 - self.i0_offset
        bstop_corrected = bstop - self.bstop_offset

        # Log transmission calculations for first runs only
        if self._log_pending:
            logger.info(f"  Calculating transmission: ratio = {transmission_ratio:.6f}")
        
        # Calculate material properties if thickness not provided
        if self.thickness is None:
            
            material_props = self.calculate_sld_mu_thickness(transmission_ratio)
            thickness = material_props['thickness_m']
            
            # Log material calculations for first runs only
            if self._log_pending:
                logger.info(f"  Material properties calculated:")
                logger.info(f"    Compound: {self.compound}")
                logger.info(f"    Energy: {self.energy_keV} keV")
                logger.info(f"    Density: {self.density_g_cm3} g/cm³")
                logger.info(f"    Absorption coefficient: {material_props['mu']:.6f} 1/m")
                logger.info(f"    SLD: {material_props['sld']:.6e}")
                logger.info(f"    Calculated thickness: {thickness:.6f} m")
        else:
            thickness = self.thickness
            if self._log_pending:
                logger.info(f"  Using provided thickness: {thickness:.6f} m")
        
        # Calculate normalization factor: trans_factor * i0_corrected
        normalization_factor = transmission_factor_raw * i0
        
        return {
            'i0_corrected': i0_corrected,
            'bstop_corrected': bstop_corrected,
            'transmission_factor': transmission_factor_raw,
            'transmission_ratio': transmission_ratio,
            'thickness': thickness,
            'normalization_factor': normalization_factor,
            'metadata_dict': metadata_dict
        }
        
    def process_saxs_file(self, raw_file_path: Path):
        """
        Process a single SAXS .raw file with corrections and 1D integration.
        
        Parameters
        ----------
        raw_file_path
            Path to the SAXS .raw file
            
        Returns
        -------
        Dict
            Dictionary containing plotting data: q, intensity, error, filename, corrections
        """
                
        # Log detailed information for first file only
        if self._log_pending:
            logger.info("")
            logger.info("=" * 60)
            logger.info("FIRST SAXS FILE PROCESSING")
            logger.info("=" * 60)
            logger.info(f"Raw file: {raw_file_path}")
            logger.info(f"Called process_saxs_file with file {raw_file_path}")
            logger.info("")
        
        # Read raw detector data
        detector_data = self.read_raw_detector_15_image(raw_file_path, 'SAXS')
        
        corrections = self.get_corrections_full(raw_file_path, "SAXS")
        
        # Log corrections for first file only
        if self._log_pending:
            self._log_pending = False
        
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
                
        add_metadata_to_dat(dat_file_path= output_path, metadata_dict=corrections['metadata_dict'])
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
        
        # Log detailed information for first file only
        if self._log_pending:
            logger.info("")
            logger.info("=" * 60)
            logger.info("FIRST WAXS FILE PROCESSING")
            logger.info("=" * 60)
            logger.info(f"Raw file: {raw_file_path}")
            logger.info(f"WAXS PONI file: {self.waxs_poni_path}")
            logger.info(f"WAXS mask file: {self.poni_path / self.mask_files['waxs'] if self.mask_files['waxs'] else 'None'}")
            logger.info(f"Detector shape: {self.detector_shapes['waxs']}")
            logger.info("")
        
        detector_data = self.read_raw_detector_15_image(raw_file_path, 'WAXS')
        corrections = self.get_corrections_full(raw_file_path, "WAXS")

        # Log corrections for first file only
        if self._log_pending:
            logger.info("=" * 60)
            self._log_pending = False

        # Create output directory and filename
        output_filename = raw_file_path.name.replace('.raw', '.dat')
        output_filename = output_filename.removeprefix("b_tassone_")
        
        output_dir = self.create_output_directory(raw_file_path, 'WAXS')
        output_path = output_dir / output_filename
        output_path = output_path.with_stem(f"{output_path.stem}_WAXS")
        
        # Perform 1D radial integration with automatic file writing
        # PyFAI will automatically generate headers with detector config and column names
        q, intensity, error = self.ai_waxs.integrate1d(
            detector_data,
            self.npt_radial,
            error_model=self.error_model,
            mask=self.waxs_mask,
            normalization_factor=corrections['normalization_factor'],
            filename=str(output_path)
        )
        add_metadata_to_dat(dat_file_path= output_path, metadata_dict=corrections['metadata_dict'])

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
    
    if (experiment.mode == "SAXS"):
        # Don't include WAXS files in just SAXS
        waxs_files = []
    else:
        waxs_pattern = str(data_directory_path / "**/WAXS/**/*.raw")
        waxs_files = [Path(f) for f in glob.glob(waxs_pattern, recursive=True)
                    if not f.endswith('.raw.pdi')]

        
    return sorted(saxs_files), sorted(waxs_files)

def full_correction_integration(config_file_path = 'config.yml', plotting=False):
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
    # Log function entry with parameters
    if (Path(config_file_path).suffix != ".yml") and (Path(config_file_path).suffix != ".yaml"):
        raise RuntimeError(f"config_file_path {config_file_path} is not a YML File. Please use a YML file.")
    print(f"full_correction_integration() called with plotting = {plotting}, config_file_path = {config_file_path}")
    print("SWAXS Data Correction and 1D Integration Demo")
    print("=" * 50)

    
    
    # Initialize experiment from config
    experiment = Experiment(config_file_path)
    logger.info(f"Calling Experiment('{config_file_path}') initialization")


    print(f"Energy: {experiment.energy_keV} keV, Density: {experiment.density_g_cm3} g/cm³")
    print()
    
    # Find all .raw files to process
    logger.info(f"Calling find_all_raw_files(experiment, {experiment.data_directory_2d})")
    saxs_files, waxs_files = find_all_raw_files(experiment, experiment.data_directory_2d)
    logger.info("")
    logger.info(f"Found {len(saxs_files)} SAXS files and {len(waxs_files)} WAXS files to process")
    logger.info("=" * 60)
    
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
            # Log function call for first run only
            if experiment._log_pending:
                logger.info(f"Calling process_saxs_file({saxs_file})")
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
            # Log function call for first run only
            if experiment._log_pending:
                logger.info(f"Calling process_waxs_file({waxs_file})")
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


def add_metadata_to_dat(dat_file_path: Path, metadata_dict: dict):
    """
    Add metadata to the end of a .dat file as commented information.
    
    Parameters
    ----------
    dat_file_path : Path
        Path to the .dat file to modify
    metadata_dict : dict
        Dictionary containing metadata with keys as column names and values as data
        
    Notes
    -----
    - Appends metadata as comments (lines starting with #) to end of .dat file
    - Writes dictionary keys and values in a structured format
    - Preserves original .dat file structure and data columns
    """
    if not dat_file_path.exists():
        raise FileNotFoundError(f"DAT file not found: {dat_file_path}")
    
    # Append metadata dictionary to the .dat file
    with open(dat_file_path, 'a') as f:
        f.write("\n# METADATA INFORMATION (YML FORMAT)\n")
        
        for key, value in metadata_dict.items():
            # Skip the empty key if it exists
            if key == "#":
                continue
            # Format the key-value pair as a comment
            f.write(f"# {key}: {value}\n")
        
def main():
    # Add command line argument
    parser = argparse.ArgumentParser(description='Process a file and generate reports')
    parser.add_argument('config_file', help='Configuration File (YML)', type = str)
    args = parser.parse_args()
    config_filename = args.config_file
    return full_correction_integration(config_file_path= config_filename, plotting=False)

if __name__ == "__main__":

    main()
